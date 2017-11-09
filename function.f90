module function
  Implicit none
  Real*8, parameter :: T=300., R=280., tau=5., dt=0.01, tf=5*tau, pi=acos(-1.)
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

  Subroutine Solve(E,E1,Tint,nbit)
    Real*8, dimension(Np), intent(inout) :: E,E1
    Real*8, dimension(0:nbit), intent(inout) :: Tint
    Integer, intent(in) :: nbit
    Integer :: i,j
    Real*8, dimension(Np) :: sig

    Do i=1,nbit
       E1=E
       Tint(i)=0
       Call sigma(sig)
       Do j=1,Np
          E(j)=(1./(1.+2*(dt/tau)))*(E1(j)+R*T*dt*(1.+sig(j)**2)/tau+&
               &2*sqrt(dt*R*T*E1(j)/tau)*sig(j))
          Tint(i)=Tint(i)+E(j)
       end Do
       Tint(i)=Tint(i)/Np
    end Do

  end Subroutine Solve

Subroutine sigma(sig)
    Real*8,dimension(Np) :: sig
    Real*8 :: u1,u2
    Integer :: i

    Do i=1,Np
       Call random_number(u1)
       Call random_number(u2)
       sig(i) = sqrt(-2*(log(u1)))*cos(2*pi*u2)
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
