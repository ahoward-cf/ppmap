        function tau2mass(kappa300,cell,d)

!--------------------------------------------------------------------
! Calculate conversion factor from optical depth at 300 microns to 
! column density in solar masses per cell.
!
! Input parameters:
!       kappa300        =       reference opacity at 300 microns [cm^2/g]
!	cell		=	sampling interval [arcsec]
!	d		=	distance [pc]
!--------------------------------------------------------------------

        implicit real(a-h,o-z)
        implicit integer(i-n)

        real, parameter :: msolar   =1.9889e33 ! solar mass [g]
        real(4) kappa300

        a = (cell*d*1.496e13)**2               ! cell area [cm^2]
        tau2mass = (a/kappa300)/msolar
        return

        end function tau2mass
