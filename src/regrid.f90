        subroutine regrid(a,ar,nx,ny,xref,yref,x0,y0,rmag)

! Regrid array a (nx by ny) onto an array, ar, of the same size, such
! that pixel (x0,y0) gets moved to (xref,yref) with a magnification
! factor of rmag.

        implicit real (a-h,o-z)
        implicit integer (i-n)

        real,   intent(in),                        &
        &       dimension(nx,ny)        ::   a      ! input array

        real,   intent(out),                       &
        &       dimension(nx,ny)        ::   ar     ! output (resampled) array

        real(8),   allocatable             ::   grid(:,:)
        real(8) x,xmin,xmax,y,ymin,ymax,z,zx,zy,zxx,zyy,zxy
        integer(4) xref,yref

        allocate( grid(0:nx, 0:ny) )
        grid = 0.
        do j = 1,ny
        do i = 1,nx
            grid(i,j) = a(i,j)
        enddo
        enddo
        do j = 1,ny
        do i = 1,nx
            x = x0 + (i-xref)/rmag     ! these are the original x and y values
            y = y0 + (j-yref)/rmag     ! sampled by i and j in the new array
            if (x >=1.d0 .and. x <= nx*1.d0 .and. &
                y >=1.d0 .and. y <= ny*1.d0) then
                call intrp2(x, 0.d0, nx*1.d0, nx, y, 0.d0, &
                    ny*1.d0, ny, grid, z,zx,zy,zxx,zyy,zxy)
            else
                z = 0.
            endif
                
            ar(i,j) = z
        enddo
        enddo

        deallocate(grid)
        return

        end subroutine regrid
