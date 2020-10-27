        subroutine resample(a,nxmod,nymod,ar,nx,ny,rmag)

! Resample array a (nxmod by nymod) onto a grid ar (nx by ny) at
! intervals of 1/rmag, mapping the centre pixel of the old array
! into the centre pixel of the new one.

        implicit real (a-h,o-z)
        implicit integer (i-n)

        real,   intent(in),                        &
        &       dimension(nxmod,nymod)  ::   a      ! input array

        real,   intent(out),                       &
        &       dimension(nx,ny)        ::   ar     ! output (resampled) array

        real(8),   allocatable             ::   grid(:,:)
        real(8) x,xmin,xmax,y,ymin,ymax,z,zx,zy,zxx,zyy,zxy

        allocate( grid(0:nxmod, 0:nymod) )
        grid = 0.

        do j = 1,nymod
        do i = 1,nxmod
            grid(i,j) = a(i,j)
        enddo
        enddo

        ixmod = nxmod/2 + 1
        iymod = nymod/2 + 1
        ix = nx/2 + 1
        iy = ny/2 + 1

        do j = 1,ny
        do i = 1,nx
            x = ixmod + (i-ix)/rmag
            y = iymod + (j-iy)/rmag
            if (x >=1.d0 .and. x <= nxmod*1.d0 .and. &
                y >=1.d0 .and. y <= nymod*1.d0) then
                call intrp2(x, 0.d0, nxmod*1.d0, nxmod, y, 0.d0, &
                    nymod*1.d0, nymod, grid, z,zx,zy,zxx,zyy,zxy)
            else
                z = 0.
            endif
                
            ar(i,j) = z
        enddo
        enddo

        deallocate(grid)
        return

        end subroutine resample
