C *************************************************************************
      subroutine writerho(outfile,array,nx,ny,nz,nbeta,ctype1,ctype2, 
     &  crpix1,crpix2,crval1,crval2,cdelt1,cdelt2,crota2,ix0,iy0,
     &  Tmin,Tmax,eta,rhoobj,Npr,Npos,dkpc,nsteps,rchisq,status)

C  Write a 4D FITS image of differential column density and associated
C  keywords.

      integer status,unit,blocksize,bitpix,naxis,naxes(4)
      integer i,j,group,fpixel,nelements
      real array(nx,ny,nz)
      real Npr,Npos
      character filename*80
      character*8 ctype1,ctype2
      character*(*) outfile
      logical simple,extend

C  The STATUS parameter must be initialized before using FITSIO.  A
C  positive value of STATUS is returned whenever a serious error occurs.
C  FITSIO uses an `inherited status' convention, which means that if a
C  subroutine is called with a positive input value of STATUS, then the
C  subroutine will exit immediately, preserving the status value. For 
C  simplicity, this program only checks the status value at the end of 
C  the program, but it is usually better practice to check the status 
C  value more frequently.

      status=0

C  Name of the FITS file to be created:
      filename=outfile

C  Delete the file if it already exists, so we can then recreate it.
C  The deletefile subroutine is listed at the end of this file.
      call deletefile(filename,status)

C  Get an unused Logical Unit Number to use to open the FITS file.
C  This routine is not required;  programmers can choose any unused
C  unit number to open the file.
      call ftgiou(unit,status)

C  Create the new empty FITS file.  The blocksize parameter is a
C  historical artifact and the value is ignored by FITSIO.
      blocksize=1
      call ftinit(unit,filename,blocksize,status)

C  Initialize parameters about the FITS image.
C  BITPIX = 16 means that the image pixels will consist of 16-bit
C  integers.  The size of the image is given by the NAXES values. 
C  The EXTEND = TRUE parameter indicates that the FITS file
C  may contain extensions following the primary array.
      simple=.true.
      bitpix=-32
      naxis=4
      naxes(1)=nx
      naxes(2)=ny
      naxes(3)=nz
      naxes(4)=nbeta
      extend=.false.

C  Write the required header keywords to the file
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
      if (status .gt. 0) then
          print *,'Error creating new FITS file.'
          if (status .eq. 105) print *,
     &        'Does a file with this name already exist?'
          stop
      endif

C  Write the array to the FITS file.
C  The last letter of the subroutine name defines the datatype of the
C  array argument ('I' = I*2, 'E' = Real*4, 'D' = Real*8).
C  The 2D array is treated as a single 1-D array with NAXIS1 * NAXIS2
C  total number of pixels.  GROUP is seldom used parameter that should
C  almost always be set = 1.
      group=1
      fpixel=1
      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
      call ftppre(unit,group,fpixel,nelements,array,status)

C  Write other optional keywords to the header
      call ftpkys(unit,'BUNIT','10^20 cm^-2',
     &  'Column density units',status)
      call ftpkys(unit,'CTYPE1',ctype1, 
     &  'name of the coordinate axis',status)
      call ftpkys(unit,'CTYPE2',ctype2, 
     &  'name of the coordinate axis',status)
      call ftpkyf(unit,'CRPIX1',crpix1,6,
     &  'coordinate system reference pixel',status)
      call ftpkyf(unit,'CRPIX2',crpix2,6,
     &  'coordinate system reference pixel',status)
      call ftpkyf(unit,'CRVAL1',crval1,6,
     &  'coordinate value at reference pixel',status)
      call ftpkyf(unit,'CRVAL2',crval2,6,
     &  'coordinate value at reference pixel',status)
      call ftpkyf(unit,'CDELT1',cdelt1,6,
     &  'coordinate increment along axis',status)
      call ftpkyf(unit,'CDELT2',cdelt2,6,
     &  'coordinate increment along axis',status)
      call ftpkyf(unit,'CROTA2',crota2,6,
     &  'coordinate system rotation angle',status)
      call ftpkyj(unit,'XORIGIN',ix0,'x origin wrt image grid',status)
      call ftpkyj(unit,'YORIGIN',iy0,'y origin wrt image grid',status)
      call ftpkyf(unit,'TMIN',Tmin,6,'Lowest temperature [K]',status)
      call ftpkyf(unit,'TMAX',Tmax,6,'Highest temperature [K]',status)
      call ftpkyf(unit,'ETA',eta,6,'Dilution',status)
      call ftpkyf(unit,'DISTANCE',dkpc,6,'Distance [kpc]',status)
      call ftpkyf(unit,'OBJCD',rhoobj,6,'Single object coldens',status)
      call ftpkye(unit,'N_PRIOR',Npr,6,'A priori no. of objects',status)
      call ftpkye(unit,'N_POST',Npos,6,'A post. no. of objects',status)
      call ftpkyj(unit,'NSTEPS',nsteps,'No. integration steps',status)
      call ftpkyf(unit,'RCHISQ',rchisq,6,'Reduced chi squared',status)

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
      return

      end
C *************************************************************************
      subroutine deletefile(filename,status)

C  A simple little routine to delete a FITS file

      integer status,unit,blocksize
      character*(*) filename

C  Simply return if status is greater than zero
      if (status .gt. 0)return

C  Get an unused Logical Unit Number to use to open the FITS file
      call ftgiou(unit,status)

C  Try to open the file, to see if it exists
      call ftopen(unit,filename,1,blocksize,status)

      if (status .eq. 0)then
C         file was opened;  so now delete it 
          call ftdelt(unit,status)
      else if (status .eq. 103)then
C         file doesn't exist, so just reset status to zero and clear errors
          status=0
          call ftcmsg
      else
C         there was some other error opening the file; delete the file anyway
          status=0
          call ftcmsg
          call ftdelt(unit,status)
      end if

C  Free the unit number for later reuse
      call ftfiou(unit, status)
      end
