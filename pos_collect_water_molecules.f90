! Original by Sandeep  D. 09-14-16 
! Program to calculate and collect H-bonds between waters
! and z-axis for palmitic acid system
!
! Adapted for POSITION_CMD files for pure water system
! by Dan   D. 06-14-17
!
! Input 
!
!  no of water molecules 
!  no of atoms in each water 
!  no of dump files 
!  fileslist  
!
! Output
!
! fort.2  ---  C_index  Avg. p2 


program orien
implicit none
integer,parameter::dp=kind(0.d0)
real(kind=dp)::dip_mag1,dip_mag2,dip_mag3,dip_mag
real(kind=dp)::dist1,dist2,xx,xy,xz,yx,yy,yz,zx,zy,zz,box(3),angle2
real(kind=dp)::angle,angle1,dipole1(3),dipole2(3),dipole3(3),dipole(3)
real(kind=dp),allocatable::coor(:,:)
!integer,allocatable::typ(:)
integer::nmol_water,natoms_water,tot_atoms,nfiles,frame_no,count_frames
integer::ifile,iframe,i,j,counter,jj
integer,allocatable::nframes(:),nskip(:),store_watid(:)

real(kind=8),allocatable :: dvec1(:,:),dvec2(:,:),dvec3(:,:),dvec4(:,:) 

real(kind=4),allocatable :: dist_factor(:),angle_factor(:)
real(kind=8),parameter :: pi=atan(1.d0)*4.d0

character(len=2)::atomtyp
character(len=2),allocatable::typ(:)

character(len=100)::file1
character(len=500),allocatable::dump_file(:)

!integer :: atomtyp
integer :: l, count

integer :: tmp_count,s,count_Ow,count1,count2,count3,count4,count5,count_collect
integer :: tmpid_Ow(1000),id_Ow(1000)


call getarg(1,file1)
open(unit=1,file=file1,action='read')  ! read from input.collect

open(unit=55555,file='WATER.MOLINDEX',action='write') 
read(1,*) nmol_water
read(1,*) natoms_water

tot_atoms=nmol_water*natoms_water
read(1,*) atomtyp
!read(1,*) atomtyp2
read(1,*) dist1
read(1,*) dist2
read(1,*) angle1
read(1,*) angle2
read(1,*) nfiles

!dist1=dist1*dist1 
!dist2=dist2*dist2 

!x=1  

allocate(store_watid(2*nmol_water))
allocate(dump_file(nfiles),nframes(nfiles),nskip(nfiles))
allocate(coor(tot_atoms,3),typ(tot_atoms))
allocate(dist_factor(2*nmol_water),angle_factor(2*nmol_water))
allocate(dvec1(nmol_water,3),dvec2(nmol_water,3),dvec3(nmol_water,3),dvec4(nmol_water,3)) 


! read dump filenames
do i=1,nfiles
   read(1,*)dump_file(i),nframes(i),nskip(i)
enddo

frame_no=0; count_frames=0
count=0 ; count1=0 ; count2=0 ; count3=0 ; count4=0 ; count5=0
counter=0
do ifile=1,nfiles
   open(unit=10+ifile,file=dump_file(ifile),action='read')

   do i=1,nskip(ifile)
     do j=1,tot_atoms+9
        read(10+ifile,*)
     enddo
   enddo

   nframes(ifile)=nframes(ifile)-nskip(ifile)

   write(*,*) 'working on file no ',ifile,' ',trim(adjustl(dump_file(ifile))),' ... '

   do iframe=frame_no+1,frame_no+nframes(ifile)
      count_frames=count_frames+1
     
!      do i=1,5
         read(10+ifile,*)   ! unneeded first line
!      enddo
 
      ! read box dimensions
      read(10+ifile,*) xx,xy,xz
      box(1)=xx       
      read(10+ifile,*) yx,yy,yz
      box(2)=yy       
      read(10+ifile,*) zx,zy,zz
      box(3)=zz       

!  read coordinates
      l=0 ; s=0 
      do i=1,tot_atoms
         read(10+ifile,*)typ(i),coor(i,:) !,coor(j,2),coor(j,3)
!         write(*,*) typ(i),coor(i,:)
      enddo
      do j=1,tot_atoms 
         if (typ(j)==atomtyp) then

             l=l+1 
             dvec1(l,:)=coor(j,:)   ! first OW
 
         endif 

         if(typ(j)==atomtyp) then

             s=s+1 
             dvec2(s,:)=coor(j,:)    ! second OW
             dvec3(s,:)=coor(j+1,:)  ! first HW of second OW
             dvec4(s,:)=coor(j+2,:)  ! second HW of second OW
         endif 

      enddo
      count_Ow=0 ; id_Ow=0 
      count_collect=0
!      write(55555,*) 'frame = ',count_frames   
    do i=1,nmol_water    ! outer Ow loop
        tmp_count=0 ; tmpid_Ow=0 
      do j=1,nmol_water    ! inner Ow loop
   
          if(i/=j) then
!            goto 2300
!          endif

          dip_mag1=0.d0;dip_mag=0.d0; angle=0.d0 
          dipole1=dvec2(j,:)-dvec1(i,:)
          dipole1(:)=dipole1(:)-box(:)*anint(dipole1(:)/box(:))
          dip_mag1=dot_product(dipole1(:),dipole1(:))
          dip_mag1=sqrt(dip_mag1) 

          if(dip_mag1<=dist1) then 
        

             dipole2=dvec2(j,:)-dvec3(j,:)          
             dipole2(:)=dipole2(:)-box(:)*anint(dipole2(:)/box(:))
             dip_mag2=sqrt(dot_product(dipole2(:),dipole2(:)))
             dipole3=dvec2(j,:)-dvec4(j,:)
             dipole3(:)=dipole3(:)-box(:)*anint(dipole3(:)/box(:))
             dip_mag3=sqrt(dot_product(dipole3(:),dipole3(:)))


             dipole(:)=dipole2(:)
             dip_mag=dip_mag2
             angle=dot_product(dipole1(:),dipole(:)/(dip_mag1*dip_mag))
             angle=dacos(angle)*180.d0/pi
             if(angle<=angle1) then 
                tmp_count=tmp_count+1 
                count_collect=count_collect+1
                tmpid_Ow(tmp_count)=j
                store_watid(count_collect)=j
                dist_factor(count_collect)=1.d0
                angle_factor(count_collect)=1.d0
!write(*,*)  j,i,dist_factor(j),count_collect
            elseif (angle<angle2 .and. angle>angle1)then
                tmp_count=tmp_count+1 
                count_collect=count_collect+1
                tmpid_Ow(tmp_count)=j
                store_watid(count_collect)=j
                dist_factor(count_collect)=1.d0
                angle_factor(count_collect)=1/(1+exp(2*(angle-32.5)))
             endif

             dipole(:)=dipole3(:)
             dip_mag=dip_mag3
             angle=dot_product(dipole1(:),dipole(:)/(dip_mag1*dip_mag))
             angle=dacos(angle)*180.d0/pi
             if(angle<=angle1) then
                tmp_count=tmp_count+1 
                count_collect=count_collect+1
                tmpid_Ow(tmp_count)=j
                store_watid(count_collect)=j
                dist_factor(count_collect)=1.d0
                angle_factor(count_collect)=1.d0
!write(*,*)  j,i,dist_factor(j),count_collect
            elseif (angle<angle2 .and. angle>angle1)then
                tmp_count=tmp_count+1 
                count_collect=count_collect+1
                tmpid_Ow(tmp_count)=j
                store_watid(count_collect)=j
                dist_factor(count_collect)=1.d0
                angle_factor(count_collect)=1/(1+exp(2*(angle-32.5)))
             endif

           
!write (*,*) dip_mag2,dip_mag3,dip_mag,dip_mag1,dist1, angle
!write (*,*) tmp_count

          elseif (dip_mag1>dist1 .and. dip_mag1<dist2) then

             dipole2=dvec2(j,:)-dvec3(j,:)          
             dipole2(:)=dipole2(:)-box(:)*anint(dipole2(:)/box(:))
             dip_mag2=sqrt(dot_product(dipole2(:),dipole2(:)))
             dipole3=dvec2(j,:)-dvec4(j,:)
             dipole3(:)=dipole3(:)-box(:)*anint(dipole3(:)/box(:))
             dip_mag3=sqrt(dot_product(dipole3(:),dipole3(:)))

!write(*,*) j,i,dip_mag1
             dipole(:)=dipole2(:)
             dip_mag=dip_mag2
             angle=dot_product(dipole1(:),dipole(:)/(dip_mag1*dip_mag))
             angle=dacos(angle)*180.d0/pi
             if(angle<=angle1) then 
                tmp_count=tmp_count+1 
                count_collect=count_collect+1
                tmpid_Ow(tmp_count)=j
                store_watid(count_collect)=j
                dist_factor(count_collect)=1/(1+exp(50*(dip_mag1-3.4)))
                angle_factor(count_collect)=1.d0

!write(*,*)  j,i,dist_factor(j),count_collect
            elseif (angle<angle2 .and. angle>angle1)then
                tmp_count=tmp_count+1 
                count_collect=count_collect+1
                tmpid_Ow(tmp_count)=j
                store_watid(count_collect)=j
                dist_factor(count_collect)=1/(1+exp(50*(dip_mag1-3.4)))
                angle_factor(count_collect)=1/(1+exp(2*(angle-32.5)))
             endif

             dipole(:)=dipole3(:)
             dip_mag=dip_mag3
             angle=dot_product(dipole1(:),dipole(:)/(dip_mag1*dip_mag))
             angle=dacos(angle)*180.d0/pi
             if(angle<angle1) then
                tmp_count=tmp_count+1 
                count_collect=count_collect+1
                tmpid_Ow(tmp_count)=j
                store_watid(count_collect)=j
                dist_factor(count_collect)=1/(1+exp(50*(dip_mag1-3.4)))
                angle_factor(count_collect)=1
!write(*,*)  j,i,dist_factor(j)
             elseif (angle<angle2 .and. angle>angle1)then
                tmp_count=tmp_count+1 
                count_collect=count_collect+1
                tmpid_Ow(tmp_count)=j
                store_watid(count_collect)=j
                dist_factor(count_collect)=1/(1+exp(50*(dip_mag1-3.4)))
                angle_factor(count_collect)=1/(1+exp(2*(angle-32.5)))
             endif

           endif
!2300  continue
         endif
        enddo       ! water loop inner

      enddo      ! water loop outer
      write(55555,*) 'frame ',count_frames,count_collect
      do jj=1,count_collect
        write(55555,*) store_watid(jj),min(dist_factor(jj),angle_factor(jj))
      enddo
   enddo        ! frame 
   frame_no=frame_no+nframes(ifile)
enddo

end program orien

