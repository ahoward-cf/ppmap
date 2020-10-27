        function planckfn(lambda,T)

! Calculate intensity [Jy/arcsec^2] at wavelength lambda [microns] of
! optically thick source of temperature T [K].
!
! Restriction:  lambda*T > 200

        real,   intent(in)    ::  lambda        ! wavelength [microns]   
        real,   intent(in)    ::  T             ! temperature [K]   

        planckfn = (977.532/lambda)**3/(exp(14401.9/(lambda*T)) - 1.)

        end function planckfn
