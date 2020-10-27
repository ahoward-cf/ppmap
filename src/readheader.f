      subroutine readheader(infile,nx,ny,crpix1,crpix2,cdelt2,status)

C  Get image size from header.

      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      integer naxes(2)
      character filename*80,comment*80,record*80
      character*(*) infile

C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

C     open the FITS file, with read-only access.  The returned BLOCKSIZE
C     parameter is obsolete and should be ignored. 
      filename = infile
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)
      if (status .eq. 222) 
     &  print *,'BITPIX must be 2nd keyword in FITS header'
      if (status .ne. 0) return

C  Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)

C  Check that it found both NAXIS1 and NAXIS2 keywords.
      if (nfound .ne. 2)then
          print *,'READHEADER failed to read the NAXISn keywords.'
          return
       end if
       nx = naxes(1)
       ny = naxes(2)

C  Read other keywords
       call ftgkye(unit,'CRPIX1',crpix1,comment,status)
       call ftgkye(unit,'CRPIX2',crpix2,comment,status)
       call ftgkye(unit,'CDELT2',cdelt2,comment,status)
       if (status .eq. 202) then
          status = 0
          call ftgkye(unit,'CD2_2',cdelt2,comment,status)
       endif

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
      return

      end
