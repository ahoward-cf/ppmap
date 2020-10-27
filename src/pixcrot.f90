        subroutine pixcrot(i,j,iref,jref,crota2,irot,jrot)

! Transform pixel coordinates into a rotated frame where the y-axis 
! corresponds to north.

        implicit real(8) (a-h, o-z)
        real(8),   parameter   :: dtor    = 0.0174533   ! degrees to radians
        real(8) i, j, iref, jref, irot, jrot

        st = sin(crota2*dtor)
        ct = cos(crota2*dtor)
        xp = (i-iref)*ct + (j-jref)*st
        yp = (j-jref)*ct - (i-iref)*st
        irot = iref + xp
        jrot = jref + yp
        return

        end subroutine pixcrot
