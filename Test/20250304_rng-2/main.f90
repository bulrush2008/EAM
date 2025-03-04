
! File: main.f90
program demo
  use rng_mod
  implicit none
  type(RNG) :: task1, task2
  integer, parameter :: N = 200000
  real :: x, sum1, sum2, sum_sq1, sum_sq2
  integer :: i

  ! Initialize with different seeds
  call task1%init(123456789_8)
  call task2%init(987654321_8)

  ! Test task1
  sum1 = 0.0; sum_sq1 = 0.0
  do i = 1, N
    !x = task1%normal() ! 高斯分布
    x = task1%uniform()
    ! debug
    !print *, "for uniform dist: x = ", x
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