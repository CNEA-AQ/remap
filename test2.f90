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

    ifile='/media/usuario/cnea_ram/WRFE/2023112400/01/wrfout_d01_2023-11-24_00_00_00' !'wrfout_d01_20XX-XX-XX_XX:00:00.nc'
!   ncks -v Times,ZS,ZNU,ZNW,XLONG,XLAT,XLONG_U,XLAT_U,XLONG_V,XLAT_V,MAPFAC_UX,MAPFAC_UY,MAPFAC_VX,MAPFAC_VY,MAPFAC_MX,MAPFAC_MY,P_HYD,P_TOP,PH,PHB,P,PB,U,V,W,T,QVAPOR,QCLOUD,QICE,QRAIN,QSNOW,QGRAUP,CLDFRA,T2,Q2,U10,V10,PSFC,RAINC,RAINNC,LANDMASK,HGT,SST,PBLH,HFX,LH,LAI,SWDOWN,GLW,LU_INDEX,ALBEDO,SEAICE,SNOW,SMOIS,SINALPHA,COSALPHA  wrfout_input subset_wrfout

    var_list=["Times","ZS","ZNU","ZNW","XLONG","XLAT","XLONG_U","XLAT_U","XLONG_V","XLAT_V","MAPFAC_UX","MAPFAC_UY","MAPFAC_VX","MAPFAC_VY","MAPFAC_MX","MAPFAC_MY","P_HYD","P_TOP","PH","PHB","P","PB","U","V","W","T","QVAPOR","QCLOUD","QICE","QRAIN","QSNOW","QGRAUP","CLDFRA","T2","Q2","U10","V10","PSFC","RAINC","RAINNC","LANDMASK","HGT","SST","PBLH","HFX","LH","LAI","SWDOWN","GLW","LU_INDEX","ALBEDO","SEAICE","SNOW","SMOIS","SINALPHA","COSALPHA"]
    typ_list=[]
    call get_grid_and_proj_from_wrf(inp_file,g,p)

    do i=1,size(var_list)
        call SCRIP_remap_field(var1,var2,g1,g2,method)
    enddo
 
   call system('echo "                                          "');  
   call system('echo "(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)"');  
   call system('echo "(!) SYSTEM CALL (!)                       "');  
   call system('echo "(!) BORRO TODOS LOS NETCDF EN DIRECTORIO! "');  
   call system('echo "(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)"');  
   call system('echo "                                          "');  
   call system('rm *.nc');  

   !***************************
   !src grids specs:
   call get_gird_and_proj_from_wrfout(ifile, g1)
   !g1%gridName='testgrid'
   !g1%proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" !'+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs'
   !g1%nx=360     !
   !g1%ny=180     !
   !g1%dx= 1.0    !
   !g1%dy= 1.0    !
   !g1%ymin=-90   !
   !g1%xmin=-180  !
   allocate(var1(g1%nx,g1%ny))
   !***************************
   !dst grid specs (silam grid):
   call get_gird_and_proj_from_wrfout(ref_file, g2)
   !g2%gridName="silamgrid"
   !g2%proj4="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
   !g2%nx=36     !g2%nx=360     !g2%nx=50      !
   !g2%ny=18     !g2%ny=180     !g2%ny=63      !
   !g2%dx=10.0   !g2%xmin=-180  !g2%xmin=-72.0 !
   !g2%dy=10.0   !g2%ymin=-90   !g2%ymin=-54.0 !
   !g2%ymin=-90  !g2%dx=1.0     !g2%dx=0.4     !
   !g2%xmin=-180 !g2%dy=1.0     !g2%dy=0.4     !
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



!=!================================================================================================
!=!get grid object from diferent files: WRF, CMAQ, GRIDDESC
!= 
!=subroutine read_GRIDDESC(griddescFile,gridName, p, g)
!=  implicit none
!=  character(200),intent(in) :: griddescFile
!=  character(*) ,intent(in)  :: gridName
!=  type(proj_type), intent(inout) :: p
!=  type(grid_type), intent(inout) :: g
!=  character(20) :: row
!=  integer :: iostat
!=   integer :: wrf_map_proj_lcc    =1 &
!=              wrf_map_proj_stere  =2 &
!=              wrf_map_proj_tmerc  =3 &
!=              wrf_map_proj_latlong=6
!=   integer :: cmaq_coordtype_latlong=1 &
!=              cmaq_coordtype_lcc    =2 &
!=              cmaq_coordtype_merc   =3 &
!=              cmaq_coordtype_stere  =4 &
!=              cmaq_coordtype_utm    =5 &
!=              cmaq_coordtype_tmerc  =8 &
!=
!=  !GRIDDESC
!=  iostat=0
!=  open(unit=2,file=griddescFile,status='old',action='read',access='sequential')
!=  do while(iostat == 0)  !loop por cada fila
!=     read(2,*,iostat=iostat) row
!=     if ( trim(row) == trim(gridname)) then
!=       g%gName=row
!=       read(2,*) p%pName,g%xmin,g%ymin,g%dx,g%dy,g%nx,g%ny !projName xorig yorig xcell ycell nrows ncols
!=       rewind(2)
!=     endif
!=     if (trim(row) == trim(p%pName)) then
!=       read(2,*) p%typ,p%alp,p%bet,p%gam,p%xcent,p%ycent
!=       iostat=1
!=     endif
!=  enddo
!=  close(2)
!=
!=  ! WRF
!=  call check(nf90_open(trim(meteo_file), nf90_read , ncid ))
!=     !grid dimensions
!=     call check (nf90_inq_dimid(ncid,'west_east'     ,  dimId   )); call check (nf90_inquire_dimension(ncid, dimId,len=g%nx ))
!=     call check (nf90_inq_dimid(ncid,'south_north'   ,  dimId   )); call check (nf90_inquire_dimension(ncid, dimId,len=g%ny ))
!=     call check (nf90_inq_dimid(ncid,'bottom_top'    ,  dimId   )); call check (nf90_inquire_dimension(ncid, dimId,len=g%nz ))
!=     call check (nf90_get_att(ncid, nf90_global, "DX", g%dx) )
!=     call check (nf90_get_att(ncid, nf90_global, "DY", g%dy) )
!=  call check(nf90_close(ncid))
!=                 CEN_LAT  = -35.00001f ;
!=                 CEN_LON  = -65.f ;
!=                 TRUELAT1 = -20.f ;
!=                 TRUELAT2 = -50.f ;
!=                 STAND_LON = -65.f ;                   !map_proj =
!=                 POLE_LAT  = 90.f ;                    !1: Lambert Conformal
!=                 POLE_LON  = 0.f ;                     !2: Polar Stereographic
!=                 MAP_PROJ  = 1 ;                       !3: Mercator
!=                                                       !6: latitude and longitude (including global)
!=  srsOut="+proj=lcc +lat_1=${truelat1} +lat_2=${truelat2} +lon_0=${stand_lon} +lat_0=${ref_lat} +a=6370000.0 +b=6370000.0 +units=m"
!=
!=  !CMAQ
!=  call check(nf90_open(trim(meteo_file), nf90_read , ncid ))
!=     !grid dimensions
!=     call check (nf90_inq_dimid(ncid,'west_east'     ,  dimId   )); call check (nf90_inquire_dimension(ncid, dimId,len=g%nx ))
!=     call check (nf90_inq_dimid(ncid,'south_north'   ,  dimId   )); call check (nf90_inquire_dimension(ncid, dimId,len=g%ny ))
!=     call check (nf90_inq_dimid(ncid,'bottom_top'    ,  dimId   )); call check (nf90_inquire_dimension(ncid, dimId,len=g%nz ))
!=     call check (nf90_get_att(ncid, nf90_global, "DX", g%dx) )
!=     call check (nf90_get_att(ncid, nf90_global, "DY", g%dy) )
!=  call check(nf90_close(ncid))
!=  !  if [ $COORDTYPE == 1 ]; then           #Geographic:
!=  !   srsOut="+proj=latlong +a=6370000.0 +b=6370000.0"
!=  !elif [ $COORDTYPE == 2 ]; then     #Lambert Conformal Conic:
!=  !   srsOut="+proj=lcc +lat_1=$P_ALP +lat_2=$P_BET +lon_0=$P_GAM +lat_0=$YCENT +a=6370000.0 +b=6370000.0 +units=m"
!=  !elif [ $COORDTYPE == 3 ]; then     #General Mercator
!=  !   srsOut="+proj=merc +lat_ts=$P_ALP +a=6370000.0 +b=6370000.0"
!=  !elif [ $COORDTYPE == 4 ]; then     #General tangent Stereografic
!=  !   srsOut="+proj=stere +lat_0=$YCENT +lon_0=$P_GAM +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
!=  !elif [ $COORDTYPE == 5 ]; then     #UTM
!=  !   echo  "proyección: 5 (Universal Transverse Mercator) no soportada en esta aplicación."; stop
!=  !elif [ $COORDTYPE == 6 ]; then     #Polar Secant Stereographic
!=  !   srsOut="+proj=stere +lat_0=$YCENT +lon_0=$P_GAM +lat_ts=lat_ts +a=6370000.0 +b=6370000.0 +k_0=1.0"
!=  !elif [ $COORDTYPE == 7 ]; then     #Equatorial Mercator
!=  !   srsOut="+proj=merc +lat_ts=$P_ALP +a=6370000.0 +b=6370000.0"
!=  !elif [ $COORDTYPE == 8 ]; then     #Transverse Mercator
!=  !   echo  "proyección: 8 (Transverse Mercator) no soportada en esta aplicación."; stop
!=  !elif [ $COORDTYPE == 9 ]; then     #Lambert Azimuthal Equal-Area
!=  !   echo  "proyección: 9 (Lambert Azimuthal Equal-Area) no soportada en esta aplicación."; stop
!=  !else
!=  !   echo  "codigo de proyección invalido. COORDTYPE"; stop
!=  !fi;
!=
!=
!=end subroutine
!=!get PROJ4 from diferent files: WRF-files, GRIDDESC, CMAQ-files
!=
!=
!=subroutine get_proj4_from_file(fileName, fileType, proj4str)
!=   
!=   implicit none
!=   character(200),intent(in)    :: fileName
!=   character( 10),intent(in)    :: fileType   ! "wrf","cmaq","griddesc"
!=   character(100),intent(inout) :: proj4str   ! proj4 string
!=   
!=   
!=   
!=   select case(fileType)
!=      case ('wrf'  )              
!=         !map_type = map_type_conserv
!=         !luse_grid_centers = .false.  !(or .true. ?) im not sure.
!=      case ('cmaq' )         
!=         !map_type = map_type_bilinear
!=         !luse_grid_centers = .true.
!=      case ('griddesc')
!=         !map_type = map_type_bicubic
!=         !luse_grid_centers = .true.
!=      case default
!=         stop '[SCRIP] Error: unknown mapping method'
!=   end select
!=   
!=
!=end subroutine
!=!================================================================================================

end program
