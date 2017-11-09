program prog
  use function
  implicit none
  Real*8, dimension(Np) :: E,E1
  Real*8, dimension(:), Allocatable :: Tint
  Integer :: nbit
  real*8 :: start,finish

  call cpu_time(start)
  nbit=int(tf/dt)
  Allocate(Tint(0:nbit))
  Call Initialize(E,Tint(0))
  E1=E
  Call Solve(E,E1,Tint,nbit)

  Call WriteTint(Tint,nbit)
  Call Writef(E)

  Deallocate(Tint)
  call cpu_time(finish)
  Print*, finish-start

end program prog
