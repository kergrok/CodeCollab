module function
  Implicit none
  Real*8, parameter :: T=300., R=280., tau=5., dt=0.01, tf=5*tau, pi2=2*acos(-1.)
  Real*8, parameter :: RTdt_tau=R*T*dt/tau, dt_tau=(1./(1.+2*(dt/tau))), cte1=sqrt(RTdt_tau)
  Integer, parameter :: Np=500000

Contains

  Subroutine Initialize(E,Tint)
    Real*8, dimension(Np), intent(out) :: E
    Real*8, intent(out) :: Tint
    Integer :: i
    Real*8 :: rand

    Tint=0.
    Do i=1,Np
       Call random_number(rand)
       E(i)= 100000+rand*100000
       Tint=Tint+E(i)
    end Do
    Tint=Tint/Np

  end Subroutine Initialize

  Subroutine Solve(E,Tint,nbit)
    Real*8, dimension(Np), intent(inout) :: E
    Real*8, dimension(0:nbit), intent(inout) :: Tint
    Integer, intent(in) :: nbit
    Integer :: i,j
    Real*8, dimension(Np) :: sig

    Do i=1,nbit
       Tint(i)=0
       Call sigma(sig)
       Do j=1,Np
          ! E(j)=(1./(1.+2*(dt/tau)))*(E(j)+R*T*dt*(1.+sig(j)**2)/tau+&
              !  &2*sqrt(dt*R*T*E(j)/tau)*sig(j))
          ! E(j)=dt_tau*(E(j)+RTdt_tau*(1.+sig(j)**2)+2*sqrt(RTdt_tau*E(j))*sig(j))
          E(j)=dt_tau*((sqrt(E(j))+cte1*sig(j))**2+RTdt_tau)
          ! Tint(i)=Tint(i)+E(j)
       end Do
       Tint(i)=sum(E)/Np
    end Do

  end Subroutine Solve

Subroutine sigma(sig)
    Real*8,dimension(Np) :: sig
    Real*8 :: u1,u2,cst
    Integer :: i

    Do i=1,Np/2
       Call random_number(u1)
       Call random_number(u2)
       cst=sqrt(-2*(log(u1)))
       sig(i) = cst*cos(pi2*u2)
       sig(i+Np/2+1)= cst*sin(pi2*u2)
    end Do
  end Subroutine sigma


  Subroutine WriteTint(Tint,nbit)
    Real*8, dimension(0:nbit), intent(inout) :: Tint
    Integer, intent(in) :: nbit
    Integer :: i

    Open(unit=10,file='Tint.dat')
    do i=0,nbit
       Write(10,*) i*dt, Tint(i)
    end do
    Close(10)

  end Subroutine WriteTint

  Subroutine Writef(E)
    Real*8, dimension(Np), intent(inout) :: E
    Integer :: i,j,k
    Real*8 :: taille
    Integer, dimension(100) :: compt

    compt=0
    taille=5*R*T/100.
    Do i=1,100
       Do k=1,Np
          If ((E(k)>=taille*(i-1)).and.(E(k)<i*taille)) then
             compt(i)=compt(i)+1
          end If
       end Do
    end Do

    Open(unit=20,file='histogram.dat')
    Do j=1,100
       Write(20,*) j*taille, compt(j)
    end Do
    Close(20)

  end Subroutine Writef

end module function
