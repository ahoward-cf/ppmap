      subroutine readimage_wcs(infile,a,nx,ny,buffer,ctype1,ctype2,
     &    crpix1,crpix2,crval1,crval2,cdelt1,cdelt2,crota2,pixel,
     &    wl,sig,status)

C  Read a FITS image and WCS keywords.

      integer status,unit,readwrite,blocksize,naxes(2),nfound
      integer group,firstpix,nbuffer,npixels,i
      real nullval,a(nx,ny),buffer(*)
      logical anynull
      character filename*80,comment*80
      character*(*) infile
      character*8 ctype1,ctype2

C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

C  Open the FITS file.
      filename=infile
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

C  Read keywords
       call ftgkys(unit,'CTYPE1',ctype1,comment,status)
       call ftgkys(unit,'CTYPE2',ctype2,comment,status)
       call ftgkye(unit,'CRPIX1',crpix1,comment,status)
       call ftgkye(unit,'CRPIX2',crpix2,comment,status)
       call ftgkye(unit,'CRVAL1',crval1,comment,status)
       call ftgkye(unit,'CRVAL2',crval2,comment,status)
       call ftgkye(unit,'CDELT1',cdelt1,comment,status)
       if (status .eq. 202) then
           status = 0
           call ftgkye(unit,'CD1_1',cdelt1,comment,status)
       endif
       call ftgkye(unit,'CDELT2',cdelt2,comment,status)
       if (status .eq. 202) then
           status = 0
           call ftgkye(unit,'CD2_2',cdelt2,comment,status)
       endif
       call ftgkye(unit,'CROTA2',crota2,comment,status)
        if (status > 0) then
            crota2 = 0.
            status = 0
        endif
       call ftgkye(unit,'PIXEL',pixel,comment,status)
       call ftgkye(unit,'WAVELEN',wl,comment,status)
       call ftgkye(unit,'SIGBACK',sig,comment,status)
        status = 0

C  Initialize variables
      npixels=nx*ny
      group=1
      firstpix=1
      nullval=-999
      nbuffer=npixels
      
      call ftgpve(unit,group,firstpix,nbuffer,nullval,
     &            buffer,anynull,status)

      k = 0
      do j = 1,ny
      do i = 1,nx
        k = k+1
        a(i,j) = buffer(k)
      enddo
      enddo

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
      return

      end
