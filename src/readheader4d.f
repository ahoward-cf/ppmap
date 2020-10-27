      subroutine readheader4d(infile,nx,ny,nz,n4,status)

C  Get 4D image size from header.

      integer status,unit,readwrite,blocksize,nkeys,nspace,hdutype,i,j
      integer naxes(4)
      character filename*80,record*80
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

C  Determine the size of the image.
      call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)

C  Check that it found NAXIS1, NAXIS2, NAXIS3, and NAXIS4 keywords.
      if (nfound .ne. 4) return
       nx = naxes(1)
       ny = naxes(2)
       nz = naxes(3)
       n4 = naxes(4)

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
      return

      end
