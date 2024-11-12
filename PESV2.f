        SUBROUTINE Main(coeffs, Rs, N, O, Ar, massN, massO, massAr,
     &    pot, dercart, terminate)

        implicit none

        integer i, terminate
        Real*8 coeffs(0:10000, 0:8), Rs(0:10000), N(3), O(3), Ar(3)
        Real*8 massN, massO, massAr, pot, dercart(3, 3), costheta, RCM
        Real*8 pot_eval
 
Cf2py   intent(in) :: coeffs, Rs, N, O, Ar, massN, massO, massAr
Cf2py   intent(out) :: pot, dercart, terminate

c       terminate flag set to 0 for off
        terminate = 0

c       Convert distances (in Angstroms to Bohr)
        do i = 1, 3
        N(i) = N(i) * 1.889725989
        O(i) = O(i) * 1.889725989
        Ar(i) = Ar(i) * 1.889725989
        end do
c       Main section

        call coords(N, O, Ar, massN, massO, RCM, costheta)
        
c       Stopping Mechanism
        if (RCM .gt. 20.0d0) then
                terminate = 1 
        endif
        
        pot = pot_eval(coeffs, RCM, costheta, Rs)
        call deriv(N, O, Ar, massN, massO, RCM, costheta, coeffs, Rs,
     &    dercart)
        end SUBROUTINE Main


        SUBROUTINE DERIV(N, O, Ar, MASSN, MASSO, RCM, costheta, COEFFM,
     &  Rs, dercart)

        implicit none

        integer i, j
        
        real*8, intent(in) :: N, O, Ar, MASSN, MASSO, costheta
        real*8, intent(in) :: COEFFM, Rs
        real*8, intent(out) :: dercart,  RCM

        double precision der, xh, xh2, pots, pot_eval, BIGR, terms, dot
        double precision smallr, dot_prod, masscn, massco, CMNO, bdlen

        dimension N(3), O(3), Ar(3), der(2), pots(2), Rs(0:10000)
        dimension COEFFM(0:10000, 0:8), BIGR(3), CMNO(3)
        dimension dercart(3, 3), terms(3), smallr(3)

        xh = 1.0d-10
        xh2 = xh*2.0d0

c        call coords(N, O, Ar, MASSN, MASSO, RCM, costheta)

        pots(1) = pot_eval(COEFFM, RCM + xh, costheta, Rs)
        pots(2) = pot_eval(COEFFM, RCM - xh, costheta, Rs)
        der(1) = (pots(1) - pots(2))/(xh2)
        pots(1) = pot_eval(COEFFM, RCM, costheta + xh, Rs)
        pots(2) = pot_eval(COEFFM, RCM, costheta - xh, Rs)
        der(2) = (pots(1) - pots(2))/(xh2)

        RCM = 0.0d0
        bdlen = 0.0d0
        do i = 1, 3
        BIGR(i) = Ar(i) - (MASSO * O(i) + MASSN * N(i))/(MASSN + MASSO)
        smallr(i) = O(i) - N(i)
        CMNO(i) = (MASSN * N(i) + MASSO * O(i))/(MASSN + MASSO)
        RCM = RCM + BIGR(i)**2
        bdlen = bdlen + smallr(i)**2
        enddo
        RCM = sqrt(RCM)
        bdlen = sqrt(bdlen)
        dot = dot_prod(BIGR, smallr)
        masscn = Massn/(massn + massO)
        massco = Masso/(massn + massO)

        do i = 1, 3

c       --------------------------------------------------------------
        terms(1) = (-BIGR(i) - masscn*smallr(i))/
     &   (RCM*bdlen)
        terms(2) = (dot/(RCM*RCM*RCM*bdlen)) *MASScn* BIGR(i)
        terms(3) = (dot/(RCM*bdlen*bdlen*bdlen))*smallr(i)
c       --------------------------------------------------------------

       
        dercart(1, i) = der(1)*(-Masscn/RCM)*BIGR(i) +
     &    der(2) * (terms(1) + terms(2) + terms(3))
c       --------------------------------------------------------------
        terms(1) = (BIGR(i) - massco*smallr(i))/(RCM*bdlen)
        terms(2) = (dot/(RCM*RCM*RCM*bdlen))*massco*BIGR(i)
        terms(3) = -(dot/(RCM*bdlen*bdlen*bdlen)) * smallr(i)       
c       --------------------------------------------------------------

        dercart(2, i) = der(1)*(-Massco/RCM)*BIGR(i) +
     &    der(2) * (terms(1) + terms(2) + terms(3))

c       --------------------------------------------------------------
        terms(1) = smallr(i)/(RCM*bdlen)
        terms(2) = -(dot/(RCM*RCM*RCM*bdlen))*BIGR(i)
        terms(3) = 0.0d0 
c       --------------------------------------------------------------

        dercart(3, i) = der(1)*(BIGR(i)/RCM) +
     &     der(2)*(terms(1) + terms(2) + terms(3))

        enddo
c write(*,*) (dercart(3, i), i = 1, 3)
        end




C       Calculates The Jocobi Coordinates        
        SUBROUTINE COORDS(N, O, Ar, MASSN, MASSO, RCM, costheta)

        IMPLICIT NONE

        INTEGER i

        Real*8, intent(in) :: N, Ar, O, MASSN, MASSO
        Real*8, intent(out) :: RCM, costheta

        DOUBLE PRECISION CM, R, MAGR, MAGAr

        DIMENSION Ar(3), N(3), O(3), CM(3), R(3)
C       ---------------------------------------------------------------
C       Setting Up Values to be altered in loop. RCM is the distance
C       between the center of mass of NO and the argon atom. Costheta
c       is the dot product of the vectors required to obtain THETA. 
        rcm = 0.0d0
        costheta = 0.0d0
C       ---------------------------------------------------------------
C       Loop over the three dimensions (x, y, z)
        do i = 1, 3
C           This line describes the center of mass of the NO molecule
            CM(i) = (N(i)*MASSN + O(i)*MASSO)/(MASSN + MASSO)
C           The bond axis direction
            R(i) = O(i) - N(i)
C           Magnitudes of bond axis and Ar vectors
            MAGR = MAGR + R(i)**2
            MAGAr = MAGAr + Ar(i)**2
C           The square of RCM        
            RCM = RCM + (Ar(i) - CM(i))**2
            costheta = costheta + R(i)*Ar(i)
        enddo

c write(*,*) RCM 
C       Need to square root quantities required
        RCM = sqrt(RCM)
        MAGR = sqrt(MAGR)
        MAGAr = sqrt(MAGAr)
C       The dot product of the vectors is divided by the product of the
C       magnitudes of the vectors to obtain the cosine of the angle
        costheta = costheta/(MAGR*MAGAr)
        return
        END

        function pot_eval(coeffs, RCM, costheta, Rs)

        implicit none

        integer i

        double precision coeffs, RCM, costheta, Rs
        double precision pot_eval, POT_EVS, ractno, ractnomin
        double precision ractnomax, ractleft, prop, plas

        dimension coeffs(0:10000, 0:8), Rs(0:10000), POT_EVS(2)
        ractno = (RCM - 2.0d0)/ 0.002d0
        ractnomin = int(ractno) + 1
        ractnomax = ractnomin + 1
        ractleft = RCM - (((ractnomin - 1) * 0.002d0) + 2.0d0)
        prop = ractleft / 0.002d0
        POT_EVS(1) = 0.0d0
        POT_EVS(2) = 0.0d0

        do i = 0, 8
        POT_EVS(1)=POT_EVS(1)+coeffs(ractnomin,i)*PLAS(i,0,costheta)
        POT_EVS(2)=POT_EVS(2)+coeffs(ractnomax,i)*PLAS(i,0,costheta)
        enddo

        pot_eval = POT_EVS(1) + (POT_EVS(2) - POT_EVS(1)) * prop
        
        return
        end

        function plas(l,m,x)
c
c Double precision Legendre polynomials
c -------------------------------------------
c calculates the modified spherical harmonic clm(x,0), where
c x lies between 1 and -1, and m lies between -l and +l

        implicit double precision(a-h,o-z)
        integer, parameter :: dp=kind(0.0d0)
        integer l,m

c m needs to be greater than or equal to zero to use this factor to
c convert associated legendre functions into spherical harmonics clm(x,0)
c
c warning - there is no check to make sure the arguments are ok!!!!!!
c

        factor = dsqrt(fact(l-m)/fact(l+m))
    
        pmm = 1.0_dp
        pll = 0.0_dp
    
        if (m.gt.0) then
                somx2 = dsqrt((1.0_dp-x)*(1.0_dp+x))
                fac=1.0_dp
        do i=1,m
            pmm = -pmm * fac * somx2
            fac = fac + 2.0_dp
        end do
         end if

        if (l.eq.m) then
                plas = factor * pmm
        else
        pmmp1 = x * dble(2*m+1) * pmm
        if (l.eq.m+1) then
                plas = factor * pmmp1
        else 
        do ll = m+2,l
        pll = ( x*dble(2*ll-1)*pmmp1-dble(ll+m-1)*pmm ) / dble(ll-m)
        pmm = pmmp1
        pmmp1 = pll
        end do
        plas = factor * pll
        end if
        end if

        return
        end function plas

        function fact(n)
ccalculates factorials to double precision
c-------------------------------------------

        implicit none

        double precision fact
        integer, intent(in) :: n
        integer nfact(0:n), i

        nfact = 0
        nfact(0) = 1
    
        do i=1,n
                nfact(i) = nfact(i-1) * i
        end do
    
        fact = dble(nfact(n))

        return
        end function fact

        SUBROUTINE propagate(N, O, Ar, deriv, massN, massO, massAr, dt,
     &  Nv, Ov, Arv)

        implicit none
        
        integer i
        double precision N, O, Ar, deriv, massN, massO, massAr, acct
        double precision acc, dt, Nv, Ov, Arv, Rf, dot_prod, accdiff

        dimension N(3), O(3), Ar(3), deriv(3, 3), acc(3, 3)
        dimension Nv(3), Ov(3), Arv(3), Rf(3), accdiff(3)
        do i = 1, 3
        acc(1, i) = -deriv(1, i)/massN
        acc(2, i) = -deriv(2, i)/massO
        acc(3, i) = -deriv(3, i)/massAr
        accdiff(i) = acc(2, i) - acc(1, i)
        Rf(i) = O(i) - N(i)
        end do
        acct = sqrt(dot_prod(accdiff, accdiff))
c write(*,*) dot_prod(Rf, accdiff), acct
        
        do i = 1, 3
        N(i) = N(i) + Nv(i) * dt + 0.5 * acc(1, i) * dt**2
        O(i) = O(i) + Ov(i) * dt + 0.5 * acc(2, i) * dt**2
        Ar(i) = Ar(i) + Arv(i) * dt + 0.5 * acc(3, i) * dt**2

        Nv(i) = Nv(i) + 0.5 * acc(1, i) * dt
        Ov(i) = Ov(i) + 0.5 * acc(2, i) * dt
        Arv(i) = Arv(i) + 0.5 * acc(3, i) * dt
        end do


        return 
        end

        function dot_prod(a, b)

        implicit none

        integer i
        double precision a, b
        double precision dot_prod

        dimension a(3), b(3)

        dot_prod = 0.0d0
        do i = 1, 3
        dot_prod = dot_prod + a(i) * b(i)
        end do
        return
        end


