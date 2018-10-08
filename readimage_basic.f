      subroutine readimage_basic(infile,a,nx,ny,buffer,status)

C  Read a FITS image.

      integer status,unit,readwrite,blocksize,naxes(2),nfound
      integer group,firstpix,nbuffer,npixels,i
      real nullval,a(nx,ny),buffer(*)
      logical anynull
      character filename*80
      character*(*) infile

C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

C  Open the FITS file.
      filename=infile
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

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
