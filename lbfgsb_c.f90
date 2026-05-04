! Copyright (c) 2013 Aubrey Barnard.  This is free software.  See
! LICENSE.txt for details.

! lbfgsb_minimize provides a simple, modern C interface to the L-BFGS-B
! FORTRAN 77 code.  It compiles together with the L-BFGS-B FORTRAN 77 code
! to create a L-BFGS-B library with a C API.
!
! This is intentionally written as an external function without a Fortran
! module wrapper.  Fortran modules generate .mod files during compilation,
! which fails when building from a read-only directory such as the Go module
! cache.  The C API is fully described in lbfgsb_c.h.

function lbfgsb_minimize( &
     ! Callbacks
     func, grad, callback_data, &
     ! Dimensionality
     dim_c, &
     ! Bounds
     bounds_control_c, lower_bounds_c, upper_bounds_c, &
     ! Parameters
     approximation_size_c, f_tolerance_c, g_tolerance_c, &
     ! Input
     initial_point_c, &
     ! Result
     min_x_c, min_f_c, min_g_c, iters_c, evals_c, &
     ! Printing, logging
     print_control_c, log_function, log_function_callback_data, &
     ! Exit status
     status_message_c, status_message_length_c) &
     result(status_c) bind(c)

  use, intrinsic :: iso_c_binding
  implicit none

  ! Signature
  type(c_funptr), intent(in), value :: func, grad, log_function
  type(c_ptr), intent(in), value :: callback_data, &
       log_function_callback_data
  integer(c_int), intent(in), value :: dim_c, approximation_size_c, &
       print_control_c, status_message_length_c
  real(c_double), intent(in), value :: f_tolerance_c, g_tolerance_c
  integer(c_int), intent(in) :: bounds_control_c(dim_c)
  real(c_double), intent(in) :: lower_bounds_c(dim_c), &
       upper_bounds_c(dim_c), initial_point_c(dim_c)
  character(c_char), intent(out) :: &
       status_message_c(status_message_length_c)
  integer(c_int), intent(out) :: iters_c, evals_c
  real(c_double), intent(out) :: min_x_c(dim_c), min_f_c, &
       min_g_c(dim_c)
  integer(c_int) :: status_c

  ! Constants (previously in module lbfgsb_entry)
  integer, parameter :: dp = kind(0d0)
  integer, parameter :: task_size = 60
  integer, parameter :: char_state_size = task_size
  integer, parameter :: bool_state_size = 4
  integer, parameter :: int_state_size = 44
  integer, parameter :: real_state_size = 29
  integer, parameter :: state_size = 14

  ! Status codes (previously enum, bind(c) in module lbfgsb_c)
  integer(c_int), parameter :: &
       LBFGSB_STATUS_SUCCESS = 0, &
       LBFGSB_STATUS_APPROXIMATE = 1, &
       LBFGSB_STATUS_WARNING = 2, &
       LBFGSB_STATUS_FAILURE = 3, &
       LBFGSB_STATUS_USAGE_ERROR = 4, &
       LBFGSB_STATUS_INTERNAL_ERROR = 5

  ! Interfaces for C callbacks and the L-BFGS-B entry point.
  ! These are defined here rather than imported from a module to avoid
  ! generating .mod files.
  interface

    function objective_function_c(dim, point, objective_function_value, &
         callback_data, status_message, status_message_length) &
         result(status) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: dim
      real(c_double), intent(in) :: point(dim)
      real(c_double), intent(out) :: objective_function_value
      type(c_ptr), intent(in), value :: callback_data
      integer(c_int), intent(in), value :: status_message_length
      character(c_char), intent(out) :: &
           status_message(status_message_length)
      integer(c_int) :: status
    end function objective_function_c

    function objective_gradient_c(dim, point, objective_function_gradient, &
         callback_data, status_message, status_message_length) &
         result(status) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      integer(c_int), intent(in), value :: dim
      real(c_double), intent(in) :: point(dim)
      real(c_double), intent(out) :: objective_function_gradient(dim)
      type(c_ptr), intent(in), value :: callback_data
      integer(c_int), intent(in), value :: status_message_length
      character(c_char), intent(out) :: &
           status_message(status_message_length)
      integer(c_int) :: status
    end function objective_gradient_c

    function log_function_c(callback_data, &
         iteration, fg_evals, fg_evals_total, step_length, &
         dim, x, f, g, &
         f_delta, f_delta_bound, g_norm, g_norm_bound) &
         result(error) bind(c)
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr), intent(in), value :: callback_data
      integer(c_int), intent(in), value :: iteration, fg_evals, &
           fg_evals_total, dim
      real(c_double), intent(in), value :: step_length, f, &
           f_delta, f_delta_bound, g_norm, g_norm_bound
      real(c_double), intent(in) :: x(dim), g(dim)
      integer(c_int) :: error
    end function log_function_c

    subroutine setulb(n, m, x, l, u, nbd, f, g, factr, pgtol, &
         wa, iwa, task, iprint, csave, lsave, isave, dsave)
      implicit none
      integer, intent(in) :: n, m, nbd(n), iprint
      double precision, intent(in) :: l(n), u(n), factr, pgtol
      character(len=60), intent(inout) :: task
      character(len=60), intent(inout) :: csave
      logical, intent(inout) :: lsave(4)
      integer, intent(inout) :: iwa(3*n), isave(44)
      double precision, intent(inout) :: x(n), f, g(n), &
           wa(2*m*n+5*n+11*m*m+8*m), dsave(29)
    end subroutine setulb

  end interface

  ! Locals (scalars before arrays)
  ! Fortran versions of arguments
  procedure(objective_function_c), pointer :: func_pointer
  procedure(objective_gradient_c), pointer :: grad_pointer
  real(dp) :: point(dim_c)
  ! Variables and memory for L-BFGS-B
  integer :: print_control
  real(dp) :: func_value, f_factor
  character(len=task_size) :: task
  character(len=char_state_size) :: char_state
  character(len=state_size) :: state
  character(len=2*task_size) :: message
  logical :: bool_state(bool_state_size)
  integer :: int_state(int_state_size), &
       working_int_memory(3 * dim_c)
  real(dp) :: grad_value(dim_c), real_state(real_state_size), &
       working_real_memory( &
       2 * approximation_size_c * dim_c + 5 * dim_c + &
       11 * approximation_size_c ** 2 + 8 * approximation_size_c &
       )

  ! Convert inputs from C types to Fortran types
  call c_f_procpointer(func, func_pointer)
  call c_f_procpointer(grad, grad_pointer)
  ! Copy initial_point_c to point because point is written to
  point = initial_point_c

  ! Start with an empty status message (fill entire string with nulls)
  status_message_c = c_null_char

  ! Translate f_tolerance to f_factor.  The convergence tolerance for
  ! the objective function is computed by the L-BFGS-B code as
  ! f_factor * epsilon(1d0) but I want to express the tolerance in
  ! terms of digits of precision, analogous to g_tolerance.
  f_factor = f_tolerance_c / epsilon(1d0)

  ! Translate print_control_c which is a zero-based version of
  ! print_control (which is possibly negative)
  print_control = print_control_c - 1

  ! Initialize the state and task
  state = 'START'
  task = state

  ! Loop to do tasks and coordinate the optimization
  do while ( &
       state == 'EVAL_FG' .or. &
       state == 'NEW_X' .or. &
       state == 'WARNING' .or. &
       state == 'START')

     ! Call L-BFGS-B code
     call setulb(dim_c, approximation_size_c, point, &
          lower_bounds_c, upper_bounds_c, bounds_control_c, &
          func_value, grad_value, &
          f_factor, g_tolerance_c, &
          working_real_memory, working_int_memory, &
          task, print_control, &
          char_state, bool_state, int_state, real_state)

     ! Interpret the returned task
     call interpret_task(task, state, message)

     ! Act on the current state
     select case (state)
     case ('EVAL_FG')
        ! Calculate function and gradient.

        ! Call objective function
        status_c = func_pointer(dim_c, point, func_value, &
             callback_data, status_message_c, status_message_length_c)
        ! Terminate optimization on any error
        if (status_c /= LBFGSB_STATUS_SUCCESS) exit

        ! Call objective function gradient
        status_c = grad_pointer(dim_c, point, grad_value, &
             callback_data, status_message_c, status_message_length_c)
        ! Terminate optimization on any error
        if (status_c /= LBFGSB_STATUS_SUCCESS) exit
     case ('WARNING')
        ! TODO handle warnings
     case ('NEW_X')
        ! Call the logging function
        status_c = call_logging_function( &
             log_function, log_function_callback_data, &
             point, func_value, grad_value, g_tolerance_c, &
             int_state, real_state, status_message_c)
        ! Terminate optimization on any error
        if (status_c /= LBFGSB_STATUS_SUCCESS) exit
     end select
  end do
  ! End optimization

  ! Return statistics
  iters_c = int_state(30)  ! Current iteration
  evals_c = int_state(34)  ! Total evaluations (each eval = [F(),G()])

  ! Analyze status and state to see how to return
  if (status_c == LBFGSB_STATUS_SUCCESS) then
     ! Objective and gradient evaluations were OK but L-BFGS-B may not
     ! be.  Regardless, take what we can from the outputs.
     min_x_c = point
     min_f_c = func_value
     min_g_c = grad_value

     ! Check for normal or problematic termination
     select case (state)
     case ('CONVERGENCE')
        status_c = LBFGSB_STATUS_SUCCESS
     case ('ABNORMAL')
        status_c = LBFGSB_STATUS_APPROXIMATE
     case ('WARNING')
        status_c = LBFGSB_STATUS_WARNING
     case ('ERROR_USAGE')
        status_c = LBFGSB_STATUS_USAGE_ERROR
     case ('ERROR_INTERNAL')
        status_c = LBFGSB_STATUS_INTERNAL_ERROR
     case default
        status_c = LBFGSB_STATUS_INTERNAL_ERROR
        message = 'Error: Unrecognized state: '//task
     end select

     ! Copy task message into status message
     call convert_f_c_string(message, status_message_c)
  else
     min_x_c = 0d0
     min_f_c = 0d0
     min_g_c = 0d0
  end if

contains

  subroutine interpret_task(task, state, message)
    implicit none
    character(len=*), intent(in) :: task
    character(len=state_size), intent(out) :: state
    character(len=*), intent(out) :: message
    integer :: cut_index

    message = ' '

    cut_index = index(task, ':') - 1
    if (cut_index == -1) cut_index = len_trim(task)

    select case (task(1:cut_index))
    case ('START')
       state = 'START'
    case ('FG', 'FG_LNSRCH', 'FG_START')
       state = 'EVAL_FG'
    case ('NEW_X')
       state = 'NEW_X'
    case ('CONVERGENCE')
       state = 'CONVERGENCE'
       message = task(14:)
    case ('ABNORMAL_TERMINATION_IN_LNSRCH')
       state = 'ABNORMAL'
       message = task
    case ('WARNING')
       state = 'WARNING'
       message = task(10:)
    case ('ERROR')
       state = 'ERROR_USAGE'
       message = task(8:)
    case default
       state = 'ERROR_INTERNAL'
       message = 'Unrecognized task: '//task
    end select
  end subroutine interpret_task

  function call_logging_function( &
       log_function_pointer_c, log_function_callback_data, &
       x, f, g, g_tolerance, &
       int_state, real_state, status_message_c) result(status_c)
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_funptr), intent(in), value :: log_function_pointer_c
    type(c_ptr), intent(in), value :: log_function_callback_data
    integer, intent(in) :: int_state(int_state_size)
    real(dp), intent(in) :: x(:), f, g(:), g_tolerance, &
         real_state(real_state_size)
    character(c_char), intent(out) :: status_message_c(:)
    integer(c_int) :: status_c
    procedure(log_function_c), pointer :: log_function_pointer
    real(dp) :: step_length, f_delta

    status_c = LBFGSB_STATUS_SUCCESS

    if (c_associated(log_function_pointer_c)) then
       call c_f_procpointer(log_function_pointer_c, &
            log_function_pointer)
       step_length = real_state(4) * real_state(14)
       f_delta = abs(real_state(2) - f) / &
            max(abs(real_state(2)), abs(f), 1d0)
       status_c = log_function_pointer( &
            log_function_callback_data, &
            int_state(30), int_state(36), int_state(34), step_length, &
            size(x), x, f, g, &
            f_delta, real_state(3), real_state(13), g_tolerance &
            )
       if (status_c /= LBFGSB_STATUS_SUCCESS) then
          call convert_f_c_string('Error: Logging function failed', &
               status_message_c)
       end if
       return
    end if
  end function call_logging_function

  subroutine convert_f_c_string(string_f, string_c)
    use, intrinsic :: iso_c_binding
    implicit none
    character(len=*), intent(in) :: string_f
    character(c_char), intent(out) :: string_c(:)
    integer :: length, i

    length = min(len_trim(string_f), size(string_c) - 1)
    forall(i = 1:length) string_c(i) = string_f(i:i)
    string_c(length+1:size(string_c)) = c_null_char
  end subroutine convert_f_c_string

end function lbfgsb_minimize