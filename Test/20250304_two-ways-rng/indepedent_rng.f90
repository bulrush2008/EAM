
module independent_rng
  implicit none
  private

  ! 定义独立随机数生成器类型
  type :: IndependentRNG
    private
    integer(8) :: state  ! 内部状态变量
  contains
    procedure :: init => rng_init      ! 初始化生成器
    procedure :: uniform => rng_uniform  ! 生成[0,1)均匀分布
    procedure :: normal => rng_normal    ! 生成高斯分布（Box-Muller）
  end type IndependentRNG

  public :: IndependentRNG

contains

  ! 初始化生成器（独立种子）
  subroutine rng_init(this, seed)
    class(IndependentRNG), intent(inout) :: this
    integer(8), intent(in) :: seed
    this%state = seed
    ! 预热生成器（避免初始相关性）
    call this%uniform()
    call this%uniform()
  end subroutine

  ! 生成[0,1)均匀分布（LCG算法）
  real function rng_uniform(this) result(u)
    class(IndependentRNG), intent(inout) :: this
    integer(8), parameter :: a = 6364136223846793005_8
    integer(8), parameter :: c = 1442695040888963407_8
    integer(8), parameter :: m = 9223372036854775808_8  ! 2^63
    this%state = mod(a * this%state + c, m)
    u = real(this%state, kind=4) / real(m, kind=4)
  end function

  ! 生成高斯分布（Box-Muller变换）
  real function rng_normal(this, mu, sigma) result(r)
    class(IndependentRNG), intent(inout) :: this
    real, intent(in) :: mu, sigma
    real :: u1, u2, mag
    real, parameter :: PI = 4.0 * atan(1.0)
    logical, save :: has_spare = .false.
    real, save :: spare

    if (has_spare) then
      r = spare * sigma + mu
      has_spare = .false.
      return
    end if

    u1 = this%uniform()
    u2 = this%uniform()
    u1 = max(u1, 1.0e-10)  ! 避免log(0)
    mag = sigma * sqrt(-2.0 * log(u1))
    r = mag * cos(2.0 * PI * u2) + mu
    spare = mag * sin(2.0 * PI * u2)
    has_spare = .true.
  end function

end module independent_rng