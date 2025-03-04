
! File: main.f90
program demo
  use rng_mod
  implicit none
  type(RNG) :: task1, task2
  integer(8), parameter :: N = 2000000
  real(IPREC) :: x, y, sum1, sum2, sum_sq1, sum_sq2
  integer :: i
  real(IPREC), allocatable, dimension(:) :: ufdata1, ufdata2

  allocate(ufdata1(1:N))
  allocate(ufdata2(1:N))

  ! Initialize with different seeds
  call task1%init(123456789_8)
  call task2%init(987654321_8)

  ! Test task1 & 2
  sum1 = 0.0; sum_sq1 = 0.0
  sum2 = 0.0; sum_sq2 = 0.0

  do i = 1, N
    x = task1%uniform()
    sum1 = sum1 + x
    sum_sq1 = sum_sq1 + x**2

    ufdata1(i) = x

    y = task2%normal()
    sum2 = sum2 + y
    sum_sq2 = sum_sq2 + y**2

    ufdata2(i) = y
  end do
  print '(a, f10.6, a, f10.6)', "Task1 - Mean:", sum1/N, "  Std:", sqrt(sum_sq1/N - (sum1/N)**2)
  print '(a, f10.6, a, f10.6)', "Task2 - Mean:", sum2/N, "  Std:", sqrt(sum_sq2/N - (sum2/N)**2)

  open(1001, file="ufdata1.csv", form="formatted")
  write(1001,*) "#uniform data"
  do i = 1, N
    write(1001, "(f15.8)") ufdata1(i)
  end do
  close(1001)

  open(1002, file="ufdata2.csv", form="formatted")
  write(1002,*) "#gaussian data"
  do i = 1, N
    write(1002, "(f15.8)") ufdata2(i)
  end do
  close(1002)

end program demo