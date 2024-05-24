program test

   use netcdf
   use PROJ  
   use SCRIP

   implicit none
   type(regular_grid_type ) :: g1,g2
   character(len=256) :: iFile, oFile
   character(len=256) :: method
   real, allocatable :: var1(:,:),var2(:,:)
   ifile='in.nc'
 
   call system('echo "                                          "');  
   call system('echo "(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)"');  
   call system('echo "(!) SYSTEM CALL (!)                       "');  
   call system('echo "(!) BORRO TODOS LOS NETCDF EN DIRECTORIO! "');  
   call system('echo "(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)"');  
   call system('echo "                                          "');  
   call system('rm *.nc');  

   !***************************
   !src grids specs:
   g1%gridName='testgrid'
   g1%proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" !'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
   g1%nx=360     !
   g1%ny=180     !
   g1%dx= 1.0    !
   g1%dy= 1.0    !
   g1%ymin=-90   !
   g1%xmin=-180  !
   allocate(var1(g1%nx,g1%ny))
   !call makeFieldTest(g1,var1)            !dummy test field
   call makeFieldTest2(g1,var1)            !sharp circle field test

   call saveArrayOnNetCDF(iFile,g1,var1)
   !***************************
   !dst grid specs (silam grid):
   g2%gridName="silamgrid"
   g2%proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
   g2%nx=36     !g2%nx=360     !g2%nx=50      !
   g2%ny=18     !g2%ny=180     !g2%ny=63      !
   g2%dx=10.0   !g2%xmin=-180  !g2%xmin=-72.0 !
   g2%dy=10.0   !g2%ymin=-90   !g2%ymin=-54.0 !
   g2%ymin=-90  !g2%dx=1.0     !g2%dx=0.4     !
   g2%xmin=-180 !g2%dy=1.0     !g2%dy=0.4     !
   allocate(var2(g2%nx,g2%ny))

   !!test Bilinear   
   method='bilinear'
   ofile='ou_'//trim(method)//'.nc'
   call SCRIP_remap_field(var1,var2,g1,g2,method)
   print*, "Bilinear Remaping succesfull!"
   call saveArrayOnNetCDF(oFile,g2,var2)

   !!test Near Neighbor
   method='distwgt'
   ofile='ou_'//trim(method)//'.nc'
   call SCRIP_remap_field(var1,var2,g1,g2,method)
   print*, "Bilinear Remaping succesfull!"
   call saveArrayOnNetCDF(oFile,g2,var2)

   !!!!test Bicubic
   method='bicubic'
   ofile='ou_'//trim(method)//'.nc'
   call SCRIP_remap_field(var1,var2,g1,g2,method)
   print*, "Bicubic Remaping succesfull!"
   call saveArrayOnNetCDF(oFile,g2,var2)

   !!!test Conservative
   method='conservative'
   ofile='ou_'//trim(method)//'.nc'
   call SCRIP_remap_field(var1,var2,g1,g2,method)
   print*, "Conserv Remaping succesfull!"
   call saveArrayOnNetCDF(oFile,g2,var2)

   !(Again to see if is faster once it has the remaping files)
   !!!!test Bicubic
   method='bicubic'
   ofile='ou_'//trim(method)//'.nc'
   call SCRIP_remap_field(var1,var2,g1,g2,method)
   print*, "Bicubic Remaping succesfull!"
   call saveArrayOnNetCDF(oFile,g2,var2)
                                           
   !!!test Conservative
   method='conservative'
   ofile='ou_'//trim(method)//'.nc'
   call SCRIP_remap_field(var1,var2,g1,g2,method)
   print*, "Conserv Remaping succesfull!"
   call saveArrayOnNetCDF(oFile,g2,var2)


   print*, "Fin de la prueba."

contains

   subroutine makeFieldTest(g,var)
       implicit none
       type(regular_grid_type ) :: g
       real                :: var(:,:)
       real(8),allocatable    :: lat(:)
       real(8),allocatable    :: lon(:)
       integer :: i,j,stat!,nx,ny
       real            :: pi=3.141593
       real            :: deg2rad=3.141593/180.0
       real            :: rad2deg=180.0/3.141593
       allocate(lat(g%ny))
       allocate(lon(g%nx))
       !allocate(var(g%nx,g%ny))
                                                                              
       lon=[ (g%xmin+i*g%dx, i=1, g%nx,1) ] !lon=c(-180:180)
       lat=[ (g%ymin+i*g%dy, i=1, g%ny,1) ] !lat=c(-90:90)
       
       ! Calculate the scalar field values using array operations
       do concurrent (i = 1:g%nx, j = 1:g%ny)
           var(i, j) = 2.0 + cos(lon(i)*deg2rad)**2 * cos(2.0*lat(j)*deg2rad)
       end do
   endsubroutine

   subroutine makeFieldTest2(g,var)
       implicit none
       type(regular_grid_type ) :: g
       real                :: var(:,:)
       real(8),allocatable    :: lat(:)
       real(8),allocatable    :: lon(:)
       integer :: i,j,stat!,nx,ny
       real :: dist
       real            :: pi=3.141593
       real            :: deg2rad=3.141593/180.0
       real            :: rad2deg=180.0/3.141593
       allocate(lat(g%ny))
       allocate(lon(g%nx))
                                                                              
       lon=[ (g%xmin+i*g%dx, i=1, g%nx,1) ] !lon=c(-180:180)
       lat=[ (g%ymin+i*g%dy, i=1, g%ny,1) ] !lat=c(-90:90)
       
       ! Calculate the scalar field values using array operations
       var=0.0
       do concurrent (i = 1:g%nx, j = 1:g%ny)
           dist=sqrt(lon(i)*lon(i)+lat(j)*lat(j))
           if ( dist < 30                           ) var(i,j)=1.0
           if ( abs(lon(i)) < 10 .and. lat(j) > -20 ) var(i,j)=0.0 
       end do
   endsubroutine


   subroutine saveArrayOnNetCDF(ncfile,g,var)
      implicit none
      character(len=*),intent(in) :: ncfile
      type(regular_grid_type ) :: g
      real                :: var(:,:)
      
      real(8),  allocatable    :: lat(:)!,y(:)
      real(8),  allocatable    :: lon(:)!,x(:)
      integer :: i,j,stat,idx
      integer :: ncid, x_dim_id,y_dim_id,var_id
      character(100) :: proj4_latlon="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

      allocate(lat(g%nx*g%ny))
      allocate(lon(g%nx*g%ny))
        
      idx=1 !cell number
      do j=0,g%ny-1
        do i=0,g%nx-1
           lon(idx)=g%xmin+0.5*g%dx+g%dx*i
           lat(idx)=g%ymin+0.5*g%dy+g%dy*j
           idx=idx+1
        enddo
      enddo                                                                    
      call proj_trans(g%proj4, proj4_latlon, lon , lat)!, g%nx*g%ny) !grid_size)

      print*,"write NETCDF"
      stat=nf90_create(ncFile, NF90_CLOBBER, ncid)
          !! Defino dimensiones
          stat=nf90_def_dim(ncid, "x" , g%nx ,   x_dim_id )
          stat=nf90_def_dim(ncid, "y" , g%ny ,   y_dim_id )
          !lat
          stat=nf90_def_var(ncid, "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
          stat=nf90_put_att(ncid, var_id, "units"          , "degrees_north"      )
          stat=nf90_put_att(ncid, var_id, "long_name"      , "latitude"           )
          !lon
          stat=nf90_def_var(ncid, "lon" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
          stat=nf90_put_att(ncid, var_id, "units"          , "degrees_east"       )
          stat=nf90_put_att(ncid, var_id, "long_name"      , "longitude"          )
          !var
          stat=nf90_def_var(ncid, "var" , NF90_FLOAT, [x_dim_id, y_dim_id], var_id )
          stat=nf90_put_att(ncid, var_id, "units"          , "g/s"                 )
          stat=nf90_put_att(ncid, var_id, "long_name"      , "var mass flux"       )
      stat=nf90_enddef(ncid)
      !Abro NetCDF outFile
      stat=nf90_open(ncFile, nf90_write, ncid)
         stat=nf90_inq_varid(ncid, 'var', var_id); stat=nf90_put_var(ncid, var_id, var(:,:)) 
         stat=nf90_inq_varid(ncid, "lat", var_id); stat=nf90_put_var(ncid, var_id, reshape(lat,[g%nx,g%ny]))!lat(:)) !
         stat=nf90_inq_varid(ncid, "lon", var_id); stat=nf90_put_var(ncid, var_id, reshape(lon,[g%nx,g%ny]))!lon(:)) !
      stat=nf90_close(ncid)
      !Cierro NetCDF outFile
   end subroutine

end program

