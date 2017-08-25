program di
 implicit none 

  character(len=300) :: file1,filewatindex,a
  integer ::  nframe,i,nmol_wat,j,frame,k,real_nmol_wat,frame1

  real(kind=8),allocatable :: dipmol(:,:),dipind(:,:),factor(:) 
  real(kind=8) :: mol(4),ind(4) 

  integer,allocatable :: molid(:)
 

 nframe=25000 
 real_nmol_wat=512 
 frame1=0

  call getarg(1,file1)
  open(unit=10,file=file1,action='read')
  
  read(10,*) real_nmol_wat
  read(10,*) nframe
  read(10,*) filewatindex

 allocate(molid(real_nmol_wat),factor(real_nmol_wat))


 allocate(dipmol(real_nmol_wat,4),dipind(real_nmol_wat,4))

  open(unit=1,file='DIPMOL_CMD',action='read')
  open(unit=2,file='DIPIND_CMD',action='read')
!  open(unit=4,file='1DH_SORT_WATER.MOLINDEX',action='read')
  open(unit=4,file=filewatindex,action='read')

  open(unit=100,file='DIPMOL_CMDtest',action='write')
  open(unit=101,file='DIPIND_CMDtest',action='write')

    nmol_wat=0 


 do i=1,nframe

    do j=1,real_nmol_wat
       read(1,*) frame,dipmol(j,:)
       read(2,*) frame,dipind(j,:)
    enddo

    read(4,*) a,frame1, nmol_wat
!write(*,*)nmol_wat
    do j=1,nmol_wat
       read(4,*) molid(j),factor(j)
    enddo       

    do j=1,real_nmol_wat

     
       do k=1,nmol_wat

          if(j==molid(k)) then 

              mol=factor(k)*dipmol(j,:)
              ind=factor(k)*dipind(j,:)

              goto 1000

          endif 
       enddo

        mol=0.d0 ; ind=0.d0 

       1000 continue 

        write(100,'(i5,4f20.12)') j,mol
        write(101,'(i5,4f20.12)') j,ind

    enddo     

 enddo

end 
