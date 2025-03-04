
! File: independent_rng.f90
module rng_mod
  implicit none
  private

  type, public :: RNG
    private
    integer(8)        :: state
    real              :: spare = 0.0
    logical           :: has_spare = .false.
  contains
    procedure :: init => seed_rng
    procedure :: normal => rng_normal
  end type RNG

  integer(8), parameter :: LCG_MULT = 6364136223846793005_8
  integer(8), parameter :: LCG_INCR = 1442695040888963407_8
  integer(8), parameter :: MODULUS = 9223372036854775807_8  ! 2^63-1

contains

  subroutine seed_rng(this, seed_val)
    class(RNG), intent(inout) :: this
    integer(8), intent(in)    :: seed_val
    this%state = mod(abs(seed_val), MODULUS)
    this%has_spare = .false.
  end subroutine seed_rng

  function rng_normal(this) result(n)
    class(RNG), intent(inout) :: this
    real :: n
    real :: u, v, s

    if (this%has_spare) then
      this%has_spare = .false.
      n = this%spare
    else
      do  ! Marsaglia polar method
        ! Generate in [-1, 1) range
        u = 2.0 * real(mod(LCG_MULT * this%state + LCG_INCR, MODULUS)) / real(MODULUS) - 1.0
        this%state = mod(LCG_MULT * this%state + LCG_INCR, MODULUS)
        v = 2.0 * real(mod(LCG_MULT * this%state + LCG_INCR, MODULUS)) / real(MODULUS) - 1.0
        this%state = mod(LCG_MULT * this%state + LCG_INCR, MODULUS)
        
        s = u*u + v*v
        if (s > 0.0 .and. s < 1.0) exit
      end do

      s = sqrt(-2.0 * log(s) / s)
      n = u * s     ! First normal value
      this%spare = v * s  ! Store second value
      this%has_spare = .true.
    end if
  end function rng_normal

end module rng_mod

! File: main.f90
program demo
  use rng_mod
  implicit none
  type(RNG) :: task1, task2
  integer, parameter :: N = 100000
  real :: x, sum1, sum2, sum_sq1, sum_sq2
  integer :: i

  ! Initialize with different seeds
  call task1%init(123456789_8)
  call task2%init(987654321_8)

  ! Test task1
  sum1 = 0.0; sum_sq1 = 0.0
  do i = 1, N
    x = task1%normal()
    sum1 = sum1 + x
    sum_sq1 = sum_sq1 + x**2
  end do
  print '(a, f10.6, a, f10.6)', "Task1 - Mean:", sum1/N, "  Std:", sqrt(sum_sq1/N - (sum1/N)**2)

  ! Test task2
  sum2 = 0.0; sum_sq2 = 0.0
  do i = 1, N
    x = task2%normal()
    sum2 = sum2 + x
    sum_sq2 = sum_sq2 + x**2
  end do
  print '(a, f10.6, a, f10.6)', "Task2 - Mean:", sum2/N, "  Std:", sqrt(sum_sq2/N - (sum2/N)**2)

end program demo