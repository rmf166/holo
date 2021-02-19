    module global

      implicit none

      integer(4),    parameter     ::                    &
          sp = kind(1.0),                                &
          dp = selected_real_kind(2*precision(1.0_sp)),  &
          qp = selected_real_kind(2*precision(1.0_dp))
      integer(4),    parameter     :: kr=qp
      integer(4)                   :: sweeps

    end module global

    program main

      implicit none

      !$ call omp_set_num_threads(4)

      call drive_fa

      contains

      subroutine drive_fa

        use global
        !$ use omp_lib, only: omp_set_num_threads

        implicit none

        integer(4)                   :: i
        integer(4)                   :: sol
        integer(4)                   :: s
        integer(4)                   :: xn
        integer(4)                   :: jmax
        integer(4),    parameter     :: n=16
        integer(4),    parameter     :: kmax=1000
        real(kind=kr)                :: x
        real(kind=kr)                :: xx
        real(kind=kr)                :: tau
        real(kind=kr)                :: c(1:4)=(/0.8000_kr,0.9000_kr,0.9900_kr,0.9999_kr/)
        real(kind=kr)                :: h(73)
        real(kind=kr)                :: rho(73,4,4)
        character(1)                 :: sw
        character(2)                 :: cols
        character(2)                 :: solopt(0:1)=(/'LD','LC'/)
        character(10)                :: caserun
        character(20)                :: datafile
        character(132)               :: datfmt

        write(cols,'(i2)') 1+4
        write(datfmt,'(a)') '(' // trim(adjustl(cols)) // '(es12.5))'

        !$ call omp_set_num_threads(n/2)

        h  =0.0_kr
        rho=0.0_kr
        do sol=0,1
          do i=1,4
            sweeps=i
            do s=1,4
              x=50.0_kr
              xx=0.82_kr
              tau=100.0_kr*(xx**46)
              do xn=1,73
                if (0.1_kr <= tau .and. tau < 1.0_kr) then
                  x=50.0_kr
                elseif (1.0_kr <= tau .and. tau < 10.0_kr) then
                  x=500.0_kr
                elseif (10.0_kr <= tau .and. tau < 500.0_kr) then
                  x=5000.0_kr
                elseif (500.0_kr <= tau .and. tau < 25000.0_kr) then
                  x=250000.0_kr
                endif
                h(xn)=tau
                jmax=x/h(xn)
                h(xn)=x/jmax
                if (mod(jmax,1) == 1) stop 'Incorrect JMAX value'
                call solve_slab_fa(sol,c(s),n,kmax,jmax,h(xn),rho(xn,s,i))
                tau=tau*(1.0_kr/xx)
              enddo
            enddo
          enddo
          do i=1,4
            write(sw,'(i1)') i
            write(caserun,'(a)') '-p1-s' // trim(adjustl(sw)) // '-' // solopt(sol)
            write(datafile,'(a)') 'numres' // trim(adjustl(caserun)) // '.dat'
            open(unit=1,file=datafile,action='write',status='unknown')
            do xn=1,73
              write(1,datfmt) h(xn),(rho(xn,s,i),s=1,4)
            enddo
            close(1)
          enddo
        enddo

      end subroutine drive_fa

      subroutine solve_slab_fa(sol,c,n,kmax,jmax,h,rho)

        use global

        implicit none

        integer(4),    intent(in)    :: sol
        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: jmax
        real(kind=kr), intent(in)    :: c
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(out)   :: rho

        integer(4)                   :: bc(2)
        real(kind=kr)                :: eps
        real(kind=kr)                :: mu(n/2)
        real(kind=kr)                :: w (n/2)
        real(kind=kr), allocatable   :: phi (:)
        real(kind=kr), allocatable   :: phi_l(:)
        real(kind=kr), allocatable   :: phi_r(:)
        real(kind=kr), allocatable   :: jnet(:)
        real(kind=kr), allocatable   :: jp(:)
        real(kind=kr), allocatable   :: jpi(:)
        real(kind=kr), allocatable   :: jm(:)
        real(kind=kr), allocatable   :: jmi(:)
        real(kind=kr)                :: q
        real(kind=kr)                :: sigt
        real(kind=kr)                :: sigs

      ! dynamic allocation of arrays

        allocate(phi(jmax))
        allocate(phi_l(jmax))
        allocate(phi_r(jmax))
        allocate(jnet(jmax+1))
        allocate(jp(jmax+1))
        allocate(jpi(jmax))
        allocate(jm(jmax+1))
        allocate(jmi(jmax))
        phi=0.0_kr
        phi_l=0.0_kr
        phi_r=0.0_kr
        jnet=0.0_kr
        jp=0.0_kr
        jpi=0.0_kr
        jm=0.0_kr
        jmi=0.0_kr

      ! build source based on options

        q=1.0_kr
        bc(1)=0
        bc(2)=1

      ! set cross sections

        sigt=1.0_kr
        sigs=c*sigt

      ! set quadrature

        call quad(n,mu,w)

      ! solve fixed-source problem

        eps=1.0e-06
        if (sol == 0) then
          call solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r,jnet,jp,jpi,jm,jmi,rho)
        elseif (sol == 1) then
          call solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r,jnet,jp,jpi,jm,jmi,rho)
        else
          write(0,'(a)') ' Incorrect solution scheme selected.'
          stop
        endif

      ! clean up arrays

        deallocate(phi)
        deallocate(phi_l)
        deallocate(phi_r)
        deallocate(jnet)
        deallocate(jp)
        deallocate(jpi)
        deallocate(jm)
        deallocate(jmi)

      end subroutine solve_slab_fa

      subroutine quad(n,mu,w)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        real(kind=kr), intent(out)   :: mu(n/2)
        real(kind=kr), intent(out)   :: w (n/2)

        integer(4)                   :: j
        integer(4),    parameter     :: nmaxp=300
        real(kind=kr)                :: xnew(nmaxp)
        real(kind=kr)                :: wnew(nmaxp)

        xnew=0.0_kr
        wnew=0.0_kr
        call gauleg(-1.0_kr,1.0_kr,xnew,wnew,n)

        do j=1,n/2
          mu(j)=xnew(j)
          w(j) =wnew(j)
        enddo

      end subroutine quad

      subroutine gauleg(x1,x2,x,w,n)
 
        use global

        implicit none

        integer(4),    intent(in)    :: n
        real(kind=kr), intent(in)    :: x1
        real(kind=kr), intent(in)    :: x2
        real(kind=kr), intent(inout) :: x(n)
        real(kind=kr), intent(inout) :: w(n)

        integer(4)                   :: i
        integer(4)                   :: j
        integer(4)                   :: m
        integer(4)                   :: kount
        integer(4),    parameter     :: nmax=300
        real(kind=kr)                :: xm
        real(kind=kr)                :: xl
        real(kind=kr)                :: p1
        real(kind=kr)                :: p2
        real(kind=kr)                :: p3
        real(kind=kr)                :: pi
        real(kind=kr)                :: pp
        real(kind=kr)                :: z
        real(kind=kr)                :: z1
        real(kind=kr)                :: xtmp(nmax)  ! full set of abscissas
        real(kind=kr)                :: wtmp(nmax)  ! full set of weights
        real(kind=kr), parameter     :: eps=1.0e-30

        pi=4.0_kr*atan(1.0_kr)
        if (n > nmax) then
          write(0,'(a,1i6)') 'Gauss-Leg. integration problem --Increase PARAMETER: NMAX to at least:',n
          stop
        endif

        m=(n+1)/2
        xm=0.5_kr*(x2+x1)
        xl=0.5_kr*(x2-x1)
        do i=1,m
          z=cos(pi*(i-0.25_kr)/(n+0.5_kr))
      1   continue
          p1=1.0_kr
          p2=0.0_kr
          do j=1,n
            p3=p2
            p2=p1
            p1=((2.0_kr*j-1.0_kr)*z*p2-(j-1.0_kr)*p3)/j
          enddo
      !   p1 is now the desired Legendre polynomial. we next compute pp, its derivative,
      !   by a standard relation involving also p2, the polynomial of one lower order.
          pp=n*(z*p1-p2)/(z*z-1.0_kr)
          z1=z
          z=z1-p1/pp
          if (abs(z-z1) > eps) go to 1
          xtmp(i)=    xm-xl*z
          xtmp(n+1-i)=xm+xl*z
      !   the (n+1-i) terms are the symmetric counterparts
          wtmp(i)=2.0_kr*xl/((1.0_kr-z*z)*pp*pp)
          wtmp(n+1-i)=wtmp(i)
        enddo

      ! (half set and assumed symmetric)
        kount=0
        do i=1,n
          if (xtmp(i) >= 0.0_kr) then
            kount=kount+1
            x(kount)=xtmp(i)   ! abscissas
            w(kount)=wtmp(i)   ! weights
          endif
        enddo

      end subroutine gauleg

      subroutine solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r,jnet,jp,jpi,jm,jmi,rho)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(inout) :: phi_l(jmax)
        real(kind=kr), intent(inout) :: phi_r(jmax)
        real(kind=kr), intent(inout) :: jnet(jmax+1)
        real(kind=kr), intent(inout) :: jp(jmax+1)
        real(kind=kr), intent(inout) :: jpi(jmax)
        real(kind=kr), intent(inout) :: jm(jmax+1)
        real(kind=kr), intent(inout) :: jmi(jmax)
        real(kind=kr), intent(inout) :: rho

        integer(4)                   :: sw
        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: kount
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: theta
        real(kind=kr)                :: norm0
        real(kind=kr)                :: norm1
        real(kind=kr)                :: psi
        real(kind=kr)                :: psil
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_out
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: alpha(:,:)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: phil(:)
        real(kind=kr), allocatable   :: philo(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: sl(:)
        real(kind=kr), allocatable   :: rhoi(:)

      ! lumping parameter (1.0_kr/3.0_kr is LD)
        theta=1.0_kr

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0_kr
        c1   =0.0_kr
        c2   =0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt*h/mu(m)
            alpha(j,m)=1.0_kr/(1.0_kr+6.0_kr/tau)
            alpha(j,m)=1.0_kr/((1.0_kr/theta-3.0_kr)*2.0_kr/tau+1.0_kr/alpha(j,m))
            c1(j,m)   =       (2.0_kr/tau+alpha(j,m)-1.0_kr)
            c2(j,m)   =1.0_kr/(2.0_kr/tau+alpha(j,m)+1.0_kr)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(sl(jmax))
        allocate(rhoi(kmax))
        phi =1.0_kr
        phil=0.0_kr
        rhoi=0.0_kr

        psi_in=0.0_kr
        psi_bc=0.0_kr
        kount=0

        allocate(phio(jmax))
        allocate(philo(jmax))
        do k=1,kmax
          phio=phi
          philo=phil
          do sw=1,sweeps
            do j=1,jmax
              s  (j)=0.5_kr*(sigs*phi  (j)+q)
              sl (j)=0.5_kr*(sigs*phil (j))
            enddo
            phi=0.0_kr
            phil=0.0_kr
            jnet=0.0_kr
            jp=0.0_kr
            jpi=0.0_kr
            jm=0.0_kr
            jmi=0.0_kr
            !$omp parallel do private(j,psi_in,psi_out,psi,psil) reduction(+:phi,phil,jp,jpi,jm,jmi)
            do m=1,n/2
              psi_in=psi_bc(m) ! left specular bc
              if (bc(1) == 0) psi_in=0.0_kr
              do j=1,jmax
                jp(j)  =jp(j)+psi_in*mu(m)*w(m)
                psi_out=c2(j,m)*(2.0_kr*(s(j)+alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
                psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr-alpha(j,m)*sl(j)/sigt
                psil   =psi_out-psi
                psi_in =psi_out
                phi(j) =phi(j)+psi*w(m)
                phil(j)=phil(j)+psil*w(m)
                jpi(j) =jpi(j)+psi*mu(m)*w(m)
              enddo
              jp(jmax+1)=jp(jmax+1)+psi_in*mu(m)*w(m)
              if (bc(2) == 0) psi_in=0.0_kr
              jm(jmax+1)=jm(jmax+1)+psi_in*mu(m)*w(m)
              do j=jmax,1,-1
                psi_out=c2(j,m)*(2.0_kr*(s(j)-alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
                psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr+alpha(j,m)*sl(j)/sigt
                psil   =psi-psi_out
                psi_in =psi_out
                phi(j) =phi(j)+psi*w(m)
                phil(j)=phil(j)+psil*w(m)
                jm(j)  =jm(j)+psi_in*mu(m)*w(m)
                jmi(j) =jmi(j)+psi*mu(m)*w(m)
              enddo
              psi_bc(m)=psi_in
            enddo
            !$omp end parallel do
            jnet=jp-jm
          enddo
          phi_l=phi-phil
          phi_r=phi+phil

          call holo_subc(sigt,sigs,h,jmax,bc,q,phi,phil,phi_l,phi_r,jp,jpi,jm,jmi)

          norm1=0.0_kr
          do j=1,jmax
            if (phi(j) > 1.0e+33 .or. phi(j) < 0.0_kr) then
              rho=10.0_kr
              return
            endif
            norm1=norm1+(phi(j)-phio(j))**2 ! +(phil(j)-philo(j))**2 is close to zero
          enddo
          norm1=sqrt(norm1)

          if (norm1 < 0.001_kr) then
            kount=kount+1
            rhoi(kount)=norm1/norm0
          endif
          norm0=norm1

          if (norm1 <= eps) exit 
          if (k > 100 .and. norm1 > 100.0_kr) exit

        enddo

        if (norm1 > eps) then
          rho=10.0_kr
        else
          rho=0.0_kr
          do j=1,kount
            rho=rho+rhoi(j)
          enddo
          rho=rho/kount
        endif

        deallocate(alpha)
        deallocate(c1)
        deallocate(c2)
        deallocate(phio)
        deallocate(phil)
        deallocate(philo)
        deallocate(s)
        deallocate(sl)
        deallocate(rhoi)

      end subroutine solve_ld

      subroutine solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,phi_l,phi_r,jnet,jp,jpi,jm,jmi,rho)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(inout) :: phi_l(jmax)
        real(kind=kr), intent(inout) :: phi_r(jmax)
        real(kind=kr), intent(inout) :: jnet(jmax+1)
        real(kind=kr), intent(inout) :: jp(jmax+1)
        real(kind=kr), intent(inout) :: jpi(jmax)
        real(kind=kr), intent(inout) :: jm(jmax+1)
        real(kind=kr), intent(inout) :: jmi(jmax)
        real(kind=kr), intent(inout) :: rho

        integer(4)                   :: sw
        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: kount
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: theta
        real(kind=kr)                :: tau3
        real(kind=kr)                :: tau5
        real(kind=kr)                :: tau7
        real(kind=kr)                :: norm0
        real(kind=kr)                :: norm1
        real(kind=kr)                :: psi
        real(kind=kr)                :: psil
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_out
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: alpha(:,:)
        real(kind=kr), allocatable   :: rbeta(:,:)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: phil(:)
        real(kind=kr), allocatable   :: philo(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: sl(:)
        real(kind=kr), allocatable   :: rhoi(:)

      ! lumping parameter (1.0_kr/3.0_kr is LD)
        theta=1.0_kr

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(rbeta(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0_kr
        rbeta=0.0_kr
        c1   =0.0_kr
        c2   =0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt*h/mu(m)
            if (tau < 0.01_kr) then
              tau3=tau *tau*tau
              tau5=tau3*tau*tau
              tau7=tau5*tau*tau
              alpha(j,m)=tau/6.0_kr-tau3/360.0_kr+tau5/15120.0_kr-tau7/604800.0_kr
            else
              alpha(j,m)=1.0_kr/tanh(tau/2.0_kr)-2.0_kr/tau
            endif
            rbeta(j,m)=1.0_kr/alpha(j,m)-6.0_kr/tau
            alpha(j,m)=1.0_kr/((1.0_kr/theta-3.0_kr)*2.0_kr/tau+1.0_kr/alpha(j,m))
            c1(j,m)=       (2.0_kr/tau+alpha(j,m)-1.0_kr)
            c2(j,m)=1.0_kr/(2.0_kr/tau+alpha(j,m)+1.0_kr)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(sl(jmax))
        allocate(rhoi(kmax))
        phi =1.0_kr
        phil=0.0_kr
        rhoi=0.0_kr

        psi_in=0.0_kr
        psi_bc=0.0_kr
        kount=0

        allocate(phio(jmax))
        allocate(philo(jmax))
        do k=1,kmax
          phio=phi
          philo=phil
          do sw=1,sweeps
            do j=1,jmax
              s  (j)=0.5_kr*(sigs*phi  (j)+q)
              sl (j)=0.5_kr*(sigs*phil (j))
            enddo
            phi=0.0_kr
            phil=0.0_kr
            jnet=0.0_kr
            jp=0.0_kr
            jpi=0.0_kr
            jm=0.0_kr
            jmi=0.0_kr
            !$omp parallel do private(j,psi_in,psi_out,psi,psil) reduction(+:phi,phil,jp,jpi,jm,jmi)
            do m=1,n/2
              psi_in=psi_bc(m) ! left specular bc
              if (bc(1) == 0) psi_in=0.0_kr
              do j=1,jmax
                jp(j)  =jp(j)+psi_in*mu(m)*w(m)
                psi_out=c2(j,m)*(2.0_kr*(s(j)+alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
                psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr-alpha(j,m)*sl(j)/sigt
                psil   =((rbeta(j,m)+1.0_kr)*psi_out+(rbeta(j,m)-1.0_kr)*psi_in)/2.0_kr-rbeta(j,m)*psi
                psi_in =psi_out
                phi(j) =phi(j)+psi*w(m)
                phil(j)=phil(j)+psil*w(m)
                jpi(j) =jpi(j)+psi*mu(m)*w(m)
              enddo
              jp(jmax+1)=jp(jmax+1)+psi_in*mu(m)*w(m)
              if (bc(2) == 0) psi_in=0.0_kr
              jm(jmax+1)=jm(jmax+1)+psi_in*mu(m)*w(m)
              do j=jmax,1,-1
                psi_out=c2(j,m)*(2.0_kr*(s(j)-alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
                psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr+alpha(j,m)*sl(j)/sigt
                psil   =-((rbeta(j,m)+1.0_kr)*psi_out+(rbeta(j,m)-1.0_kr)*psi_in)/2.0_kr+rbeta(j,m)*psi
                psi_in =psi_out
                phi(j) =phi(j)+psi*w(m)
                phil(j)=phil(j)+psil*w(m)
                jm(j)  =jm(j)+psi_in*mu(m)*w(m)
                jmi(j) =jmi(j)+psi*mu(m)*w(m)
              enddo
              psi_bc(m)=psi_in
            enddo
            !$omp end parallel do
            jnet=jp-jm
          enddo
          phi_l=phi-phil
          phi_r=phi+phil

          call holo_subc(sigt,sigs,h,jmax,bc,q,phi,phil,phi_l,phi_r,jp,jpi,jm,jmi)

          norm1=0.0_kr
          do j=1,jmax
            if (phi(j) > 1.0e+33 .or. phi(j) < 0.0_kr) then
              rho=10.0_kr
              return
            endif
            norm1=norm1+(phi(j)-phio(j))**2 ! +(phil(j)-philo(j))**2 is close to zero
          enddo
          norm1=sqrt(norm1)

          if (norm1 < 0.001_kr) then
            kount=kount+1
            rhoi(kount)=norm1/norm0
          endif
          norm0=norm1

          if (norm1 <= eps) exit 
          if (k > 100 .and. norm1 > 100.0_kr) exit

        enddo

        if (norm1 > eps) then
          rho=10.0_kr
        else
          rho=0.0_kr
          do j=1,kount
            rho=rho+rhoi(j)
          enddo
          rho=rho/kount
        endif

        deallocate(alpha)
        deallocate(rbeta)
        deallocate(c1)
        deallocate(c2)
        deallocate(phio)
        deallocate(phil)
        deallocate(philo)
        deallocate(s)
        deallocate(sl)
        deallocate(rhoi)

      end subroutine solve_lc

      subroutine holo_subc(sigt,sigs,h,jmax,bc,q,phi,phil,phi_l,phi_r,jp,jpi,jm,jmi)

        use global

        implicit none

        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(inout) :: phil(jmax)
        real(kind=kr), intent(in)    :: phi_l(jmax)
        real(kind=kr), intent(in)    :: phi_r(jmax)
        real(kind=kr), intent(in)    :: jp(jmax+1)
        real(kind=kr), intent(in)    :: jpi(jmax)
        real(kind=kr), intent(in)    :: jm(jmax+1)
        real(kind=kr), intent(in)    :: jmi(jmax)

        integer(4)                   :: j
        integer(4)                   :: n
        real(kind=kr)                :: dc
        real(kind=kr)                :: db
        real(kind=kr)                :: siga
        real(kind=kr)                :: s(2*jmax)
        real(kind=kr)                :: ccp(2*jmax+1)
        real(kind=kr)                :: ccm(2*jmax)
        real(kind=kr)                :: a(2*jmax)
        real(kind=kr)                :: b(2*jmax)
        real(kind=kr)                :: c(2*jmax)
        real(kind=kr)                :: x(2*jmax)

        if (mod(jmax,1) /= 0) stop ' Fine mesh does not align with HOLO_SUBC.'
        if (bc(1) /= 0) stop ' Only vacuum BC supported on left edge.'
        if (bc(2) /= 1) stop ' Only reflective BC supported on right edge.'

        s=0.0_kr
        ccp=0.0_kr
        ccm=0.0_kr
        a=0.0_kr
        b=0.0_kr
        c=0.0_kr
        x=0.0_kr

        dc=1.0_kr/(3.0_kr*(0.5_kr*h)*sigt)
        db=dc/(2.0_kr+2.0_kr*dc)
        do j=1,jmax
          if (j > 1) then
            ccp(2*j-1)=(2.0_kr*jp(j)+dc*(phi_l(j)-phi_r(j-1)))/(2.0_kr*phi_r(j-1))
            ccm(2*j-1)=(2.0_kr*jm(j)-dc*(phi_l(j)-phi_r(j-1)))/(2.0_kr*phi_l(j))
          else
            ccm(2*j-1)=(2.0_kr*jm(j)-db*(phi_l(j)           ))/(2.0_kr*phi_l(j))
          endif
          ccp(2*j)=(2.0_kr*jpi(j)+dc*(phi_r(j)-phi_l(j)))/(2.0_kr*phi_l(j))
          ccm(2*j)=(2.0_kr*jmi(j)-dc*(phi_r(j)-phi_l(j)))/(2.0_kr*phi_r(j))
        enddo
        siga=sigt-sigs
        b(1)=dc+ccp(2)+0.5_kr*db+ccm(1)+(0.5_kr*h)*siga
        c(1)=-(dc+ccm(2))
        do n=2,2*jmax-1
          a(n)=-(dc+ccp(n))
          b(n)=2.0_kr*dc+ccm(n)+ccp(n+1)+(0.5_kr*h)*siga
          c(n)=-(dc+ccm(n+1))
        enddo
        n=2*jmax
        a(n)=-(dc+ccp(n))
        b(n)=dc+ccm(n)+(0.5_kr*h)*siga

        s=(0.5_kr*h)*q
        x=s
        call tdma(2*jmax,a,b,c,x)

        do j=1,jmax
          phi (j)=0.5_kr*(x(2*j)+x(2*j-1))
          phil(j)=0.5_kr*(x(2*j)-x(2*j-1))
        enddo

      end subroutine holo_subc

      subroutine tdma(n,a,b,c,d)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        real(kind=kr), intent(in)    :: a(n)
        real(kind=kr), intent(inout) :: b(n)
        real(kind=kr), intent(in)    :: c(n)
        real(kind=kr), intent(inout) :: d(n)

        integer(4)                   :: k
        integer(4)                   :: km1
        real(kind=kr)                :: xm

        ! simple solution for n = 1

        if (n == 1) then
          d(1) = d(1)/b(1)
          return
        endif

        ! forward substitution

        do k=2,n
          km1 =k-1
          xm  =a(k)/b(km1)
          b(k)=b(k)-xm*c(km1)
          d(k)=d(k)-xm*d(km1)
        enddo

        ! back substitition

        d(n)=d(n)/b(n)
        do k=n-1,1,-1
          d(k)=(d(k)-c(k)*d(k+1))/b(k)
        enddo

      end subroutine tdma

    end program main
