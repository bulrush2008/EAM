
program main
  use independent_rng
  implicit none
  type(IndependentRNG) :: task1, task2
  integer, parameter :: N = 100000
  real :: data1(N), data2(N)
  integer :: i

  ! 初始化独立生成器（不同种子）
  call task1%init(seed=42_8)
  call task2%init(seed=12345_8)

  ! 交替生成随机数
  do i = 1, N
    data1(i) = task1%normal(mu=0.0, sigma=1.0)  ! 任务1：标准正态分布
    data2(i) = task2%normal(mu=5.0, sigma=2.0)  ! 任务2：μ=5, σ=2
  end do

  ! 统计验证
  print *, "Task1 - Mean:", sum(data1)/N, "Std:", sqrt(sum((data1 - sum(data1)/N)**2)/(N-1))
  print *, "Task2 - Mean:", sum(data2)/N, "Std:", sqrt(sum((data2 - sum(data2)/N)**2)/(N-1))

end program main