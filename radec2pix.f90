        subroutine radec2pix(ra,dec,i,j,ra0,dec0,iref,jref,pixel,inverse)

! Transformation for tangent-plane projection.
!
! Input parameters:
!	ra,dec		= equatorial position [deg]
!	i,j		= pixel location
!	ra0,dec0	= reference position [deg]
!	iref,jref	= reference pixel
!	pixel		= pixel size [deg]
!	inverse		= if .true., take the pixel location and calculate the
!			  ra and dec.

        implicit real(8) (a-h,o-z)

        real(8),   parameter   :: dtor    = 0.0174533   ! degrees to radians
        real(8) i, j, iref, jref, line
        logical inverse

        scale = 1.d0/pixel
        alpha0 = ra0*dtor
        delta0 = dec0*dtor

        if (inverse) then 
	    sample = i - iref
	    line = jref - j
	    x = sample*dtor/scale
	    y = line*dtor/scale
	    D = atan(sqrt(x**2 + y**2))
	    B = atan2(-x,y)
	    xx = sin(delta0)*sin(D)*cos(B) + cos(delta0)*cos(D)
	    yy = sin(D)*sin(B)
	    alpha = alpha0 + atan2(yy,xx)
	    delta = asin(sin(delta0)*cos(D) - cos(delta0)*sin(D)*cos(B))
	    ra = alpha/dtor
	    dec = delta/dtor
        else
	    alpha = ra*dtor
	    delta = dec*dtor
	    A = cos(delta)*cos(alpha - alpha0)
	    F = (scale/dtor)/(sin(delta0)*sin(delta) + A*cos(delta0))
	    line = -F*(cos(delta0)*sin(delta) - A*sin(delta0))
	    sample = -F*cos(delta)*sin(alpha - alpha0)
	    i = iref + sample
	    j = jref - line
        endif

        return

        end
