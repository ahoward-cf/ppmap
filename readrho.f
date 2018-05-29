      subroutine readrho(infile,a,nx,ny,nz,n4,buffer,ix0,iy0,eta,rhoobj,
     &  Nprior,Nexp,nsteps,rchisq,status)

C  Read a 4D FITS image of differential column density and associated
C  keywords.

      integer status,unit,readwrite,blocksize,naxes(4),nfound
      integer group,firstpix,nbuffer,npixels,i
      real nullval,a(nx,ny,nz,n4),buffer(*),Nprior,Nexp
      logical anynull
      character filename*80,comment*80
      character*(*) infile

C  The STATUS parameter must always be initialized.
      status=0

C  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

C  Open the FITS file.
      filename=infile
      readwrite=0
      call ftopen(unit,filename,readwrite,blocksize,status)

C  Read keywords.
      call ftgkyj(unit,'XORIGIN',ix0,comment,status)
      call ftgkyj(unit,'YORIGIN',iy0,comment,status)
      call ftgkye(unit,'ETA',eta,comment,status)
      call ftgkye(unit,'OBJCD',rhoobj,comment,status)
      call ftgkye(unit,'N_PRIOR',Nprior,comment,status)
      call ftgkye(unit,'N_POST',Nexp,comment,status)
      call ftgkyj(unit,'NSTEPS',nsteps,comment,status)
      call ftgkye(unit,'RCHISQ',rchisq,comment,status)
      status = 0

C  Initialize variables
      npixels=nx*ny*nz*n4
      group=1
      firstpix=1
      nullval=-999
      nbuffer=npixels
      
      call ftgpve(unit,group,firstpix,nbuffer,nullval,
     &            buffer,anynull,status)

      icount = 0
      do l = 1,n4
      do k = 1,nz
      do j = 1,ny
      do i = 1,nx
        icount = icount+1
        a(i,j,k,l) = buffer(icount)
      enddo
      enddo
      enddo
      enddo

C  The FITS file must always be closed before exiting the program. 
C  Any unit numbers allocated with FTGIOU must be freed with FTFIOU.
      call ftclos(unit, status)
      call ftfiou(unit, status)
      return

      end
