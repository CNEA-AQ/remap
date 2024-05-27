program test

   use netcdf
   use PROJ
   use SCRIP

  implicit none
  type(regular_grid_type ) :: g1,g2        !src dst
  character(len=256)       :: method
  integer :: iostat

  character(100), parameter :: proj_latlon = "+proj=latlong +a=6370000.0 +b=6370000.0"

  !dst:
  character(len=256) :: iFile, oFile !src dst
  integer            :: nx,ny,nz,time, soil_levels
  real               :: dx,dy,clon,clat
real(8) :: xcent, ycent
  character(100)     :: proj
                                                                          
  namelist/grid_specs/iFile,oFile, nx,ny,nz,soil_levels, dx,dy,clon,clat, proj

  !Leo namelist:
  read(*,nml=grid_specs, iostat=iostat)
  if( iostat /= 0 ) stop 'error reading namelist'

   !call system('rm *.nc');  

   !***************************
   !src grids specs:
   print*, "src grids specs:"
   call get_wrf_grid_and_proj(trim(iFile), g1, time, soil_levels)
   !***************************
   !dst grid specs (from namelist):
   print*, "dst grid specs (from namelist):"
   g2%gridName="Dst grid"
   g2%proj4=proj
   g2%nx=nx;   g2%ny=ny;   g2%nz=nz;
   g2%dx=dx;   g2%dy=dy
   xcent=clon; ycent=clat

   call proj_trans(proj_latlon, trim(proj), xcent, ycent )
   g2%xmin=xcent - g2%dx * 0.5 * g2%nx
   g2%ymin=ycent - g2%dy * 0.5 * g2%ny
   g2%xmax=xcent + g2%dx * 0.5 * g2%nx
   g2%ymax=ycent + g2%dy * 0.5 * g2%ny
   
!debug:
call print_grid_values(g1)
call print_grid_values(g2)


!   !allocate(var1(g1%nx,g1%ny))
!   !allocate(var2(g2%nx,g2%ny))
!
!   !***************************
!   !Init outFile
!   call create_wrf_out_file(oFile, dstFile)
!    
!
!   !***************************
!   !Remap variables:
!
!   call check(nf90_open(trim(inp_file), nf90_read, ncid ))
!
!   do i=1,size(var_list)
!
!      gridDims = grid_dims(i)   ! "2D","3D" (+ time dim)
!      gridType = grid_type(i)   ! (M,U,V,W,S)_GRID
!       varType =  var_type(i)   ! float, int
!       varName =  var_list(i)   !
!
!      g1=wrf%grd                !current grid
!
!      if ( grid_type(i) == "M-GRID" )  g1%nx=g1%nx+1; g1%xmin=g1%xmin-g1%dx; g1%ymin=g1%ymin-g1%dy; 
!      if ( grid_type(i) == "U-GRID" )  g1%nx=g1%nx+1; g1%xmin=g1%xmin-g1%dx; g1%ymin=g1%ymin-g1%dy
!      if ( grid_type(i) == "V-GRID" )  g1%ny=g1%nx+1; g1%xmin=g1%xmin-g1%dx; g1%ymin=g1%ymin-g1%dy
!      if ( grid_type(i) == "W-GRID" )  g1%nz=g1%nx+1; g1%xmin=g1%xmin-g1%dx; g1%ymin=g1%ymin-g1%dy
!
!      allocate(tmp(nx,ny,nz,Times))
!
!      call check( nf90_inq_varid(ncid,trim(varName), varid ))
!
!      do t=1,wrf%Times
!
!           select case(type_list(i))
!             case ("2D")
!                  call check(nf90_get_var (ncid, varid, tmp(:,:,t), start=(/1,1,t/),count=(/nx,ny,1/)) ) 
!                  call SCRIP_remap_field(var1,var2,g1,g2,method)
!
!             case ("3D")
!                  do z=1,grid%nz
!                     call check(nf90_get_var (ncid, varid, tmp(:,:,z,t), start=(/1,1,z,t/),count=(/nx,ny,1,1/)) ) 
!                     call SCRIP_remap_field(var1,var2,g1,g2,method)
!                  enddo
!             case default
!
!      enddo
!
!      !Save variable in outFile:
!      call save_variable_in_netcdf(oFile, var2, g2)
!
!      deallocate(tmp)
!   enddo
!
!   call check(nf90_close(ncid))
!
!   print*, "Fin de la prueba."
!
contains
 subroutine check(status)
   integer, intent(in) :: status
   if (status /= nf90_noerr) then
     write(*,*) nf90_strerror(status)
     stop 'netcdf error'
   end if
 end subroutine check

subroutine print_grid_values(g)
 implicit none
 type(regular_grid_type) , intent(in) :: g
real (8) :: minlat,maxlat,minlon,maxlon
 !debug:
 print '(A10)', g%gridName
 print '("proj:", A110)', g%proj4
 print '("nx=",I4," ny=",I4," nz=", I4," soil_layers=", I4)',g%nx, g%ny, g%nz,4
 print '("dx=",F9.3," dy=",F9.3)' ,g%dx, g%dy
 print '("xmin=",F12.3," ymin=",F12.3,/," xmax=",F12.3," ymax=",F12.3)', g%xmin,g%ymin,g%xmax, g%ymax

minlat=g%ymin;maxlat=g%ymax
minlon=g%xmin;maxlon=g%xmax
call proj_trans(g%proj4, proj_latlon,minlon,minlat ) 
call proj_trans(g%proj4, proj_latlon,maxlon, maxlat) 
 print '("lon-min=",F12.3," lat-min=",F12.3,/," lon-max=",F12.3," lat-max=",F12.3)', minlon,minlat,maxlon, maxlat

end subroutine

   subroutine get_wrf_grid_and_proj(wrf_file,g,Time, soil_layers_stag)
       implicit none

       character(len=*) :: wrf_file
       type(regular_grid_type ), intent(inout) :: g
       integer, intent(inout)                  :: Time, soil_layers_stag
       !proj
       real    :: cen_lat, cen_lon, truelat1,truelat2,stand_lon
       real    :: a=6370000.0,b=6370000.0      !axis and semi-axis radius of wrf-spheroid
       character(len=256) :: proj_str
       integer :: proj_id
       real(8) :: xcent, ycent
integer :: ncid,dimid
       !grid
       !integer :: south_north,west_east,bottom_top,time
       !integer :: dx,dy 

       call check(nf90_open(trim(wrf_file), nf90_write, ncid ))
            !Projection parameters
            call check(nf90_get_att(ncid, NF90_GLOBAL, "MAP_PROJ",  proj_id  )) !1: lcc  2: stere 3: merc  6: latlon
            call check(nf90_get_att(ncid, NF90_GLOBAL, "CEN_LAT" ,  CEN_LAT  ))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "CEN_LON" ,  CEN_LON  ))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1",  TRUELAT1 ))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2",  TRUELAT2 ))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", STAND_LON))

            !Grid dimensions:
            !DateStrLen          = 19;
            call check(nf90_inq_dimid(ncid,"Time"            ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=Time             ))
            call check(nf90_inq_dimid(ncid,"west_east"       ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=g%nx             ))
            call check(nf90_inq_dimid(ncid,"south_north"     ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=g%ny             ))
            call check(nf90_inq_dimid(ncid,"bottom_top"      ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=g%nz             ))
            !call check(nf90_inq_dimid(ncid,"west_east_stag"  ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%west_east_stag  ))
            !call check(nf90_inq_dimid(ncid,"south_north_stag",dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%south_north_stag))
            !call check(nf90_inq_dimid(ncid,"bottom_top_stag" ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%bottom_top_stag ))
            call check(nf90_inq_dimid(ncid,"soil_layers_stag",dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=soil_layers_stag ))
       
            call check(nf90_get_att(ncid, NF90_GLOBAL, "DX", g%dx))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "DY", g%dy))
       call check(nf90_close(ncid))
     
       !set Projection       
       select case (proj_id)
         case (1) !lcc
            write(proj_str,'("+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=",F8.4," +lat_2=",F8.4," +lat_0=",F8.4," +lon_0=",F8.4)') truelat1,truelat2,cen_lat,stand_lon
         case (2) !stere
            write(proj_str,'("+proj=stere +a=",F12.1," +b=",F12.1," +lat_0=",F9.4," +lon_0=",F9.4," +lat_ts=",F9.4)') a,b,truelat1,truelat2,cen_lat,stand_lon
         case (3) !merc
            write(proj_str,'("+proj=merc +a=",F12.1," +b=",F12.1," +lon_0=",F9.4," +lat_ts=",F9.4)') a,b,truelat1,truelat2,stand_lon,cen_lat
         case (6)
           write(proj_str,'("+proj=latlong +a=",F12.1," b=",F12.1)') a,b
       end select

       g%proj4=proj_str

       !set Grid parameters:
       xcent=CEN_LON; ycent=CEN_LAT
       call proj_trans(proj_latlon, proj_str, xcent, ycent )
       !M-GRID
       g%xmin=xcent - g%dx *0.5* g%nx
       g%ymin=ycent - g%dy *0.5* g%ny
       g%xmax=xcent + g%dx *0.5* g%nx
       g%ymax=ycent + g%dy *0.5* g%ny
       
   end subroutine


!   subroutine create_wrf_outfile(fileName, ncfile,g,var)
!      implicit none
!      character(len=*),intent(in) :: ncfile
!      type(regular_grid_type ) :: g
!      real                :: var(:,:)
!      
!      real(8),  allocatable    :: lat(:)!,y(:)
!      real(8),  allocatable    :: lon(:)!,x(:)
!      integer :: i,j,stat,idx
!      integer :: ncid, x_dim_id,y_dim_id,var_id
!      character(100) :: proj4_latlon="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
!
!      allocate(lat(g%nx*g%ny))
!      allocate(lon(g%nx*g%ny))
!        
!      idx=1 !cell number
!      do j=0,g%ny-1
!        do i=0,g%nx-1
!           lon(idx)=g%xmin+0.5*g%dx+g%dx*i
!           lat(idx)=g%ymin+0.5*g%dy+g%dy*j
!           idx=idx+1
!        enddo
!      enddo                                                                    
!      call proj_trans(g%proj4, proj4_latlon, lon , lat)!, g%nx*g%ny) !grid_size)
!
!      print*,"write NETCDF"
!      stat=nf90_create(ncFile, NF90_CLOBBER, ncid)
!          !! Defino dimensiones
!          stat=nf90_def_dim(ncid, "x" , g%nx ,   x_dim_id )
!          stat=nf90_def_dim(ncid, "y" , g%ny ,   y_dim_id )
!          stat=nf90_def_dim(ncid, "DateStrLen"      , 19          ,dim_id_dateStrLen      )
!          stat=nf90_def_dim(ncid, "Time"            , g%Times     ,dim_id_Time            )
!          stat=nf90_def_dim(ncid, "south_north"     , g%ny        ,dim_id_south_north     )
!          stat=nf90_def_dim(ncid, "west_east"       , g%nx        ,dim_id_west_east       )
!          stat=nf90_def_dim(ncid, "bottom_top"      , g%nz        ,dim_id_bottom_top      )
!          stat=nf90_def_dim(ncid, "west_east_stag"  , g%nx+1      ,dim_id_west_east_stag  )
!          stat=nf90_def_dim(ncid, "south_north_stag", g%ny+1      ,dim_id_south_north_stag)
!          stat=nf90_def_dim(ncid, "bottom_top_stag" , g%nz+1      ,dim_id_bottom_top_stag )
!          stat=nf90_def_dim(ncid, "soil_layers_stag", g%nSoil     ,dim_id_soil_layers_stag)
!
!          !
!          stat=nf90_def_var(ncid, "Times"  , "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "ZS"     , "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "ZNU"    , "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "ZNW"    , "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "XLONG"  , "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "XLAT"   , "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "XLONG_U", "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "XLAT_U" , "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "XLONG_V", "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_def_var(ncid, "XLAT_V" , "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!
!          !!Coordinates:
!          !allocate( Times(Time)) ;
!          !allocate(    XLAT(south_north, west_east     )) ;!Time,
!          !allocate(  XLAT_U(south_north, west_east_stag)) ;!Time,
!          !allocate(  XLAT_V(south_north_stag, west_east)) ;!Time,
!          !allocate(   XLONG(south_north, west_east     )) ;!Time,
!          !allocate( XLONG_U(south_north, west_east_stag)) ;!Time,
!          !allocate( XLONG_V(south_north_stag, west_east)) ;!Time,
!          !allocate(     ZNU(bottom_top                 )) ;!Time,
!          !allocate(     ZNW(bottom_top_stag            )) ;!Time,
!          !allocate(      ZS(soil_layers_stag           )) ;!Time,
!
!          !call check( nf90_inq_varid(ncid,"Times"  , varid )); call check( nf90_get_var(ncid, varid , Times    ))  !"Times"
!          !call check( nf90_inq_varid(ncid,"ZS"     , varid )); call check( nf90_get_var(ncid, varid , ZS       ))  !"ZS"
!          !call check( nf90_inq_varid(ncid,"ZNU"    , varid )); call check( nf90_get_var(ncid, varid , ZNU      ))  !"ZNU"
!          !call check( nf90_inq_varid(ncid,"ZNW"    , varid )); call check( nf90_get_var(ncid, varid , ZNW      ))  !"ZNW"
!          !call check( nf90_inq_varid(ncid,"XLONG"  , varid )); call check( nf90_get_var(ncid, varid , XLONG    ))  !"XLONG"
!          !call check( nf90_inq_varid(ncid,"XLAT"   , varid )); call check( nf90_get_var(ncid, varid , XLAT     ))  !"XLAT"
!          !call check( nf90_inq_varid(ncid,"XLONG_U", varid )); call check( nf90_get_var(ncid, varid , XLONG_U  ))  !"XLONG_U"
!          !call check( nf90_inq_varid(ncid,"XLAT_U" , varid )); call check( nf90_get_var(ncid, varid , XLAT_U   ))  !"XLAT_U"
!          !call check( nf90_inq_varid(ncid,"XLONG_V", varid )); call check( nf90_get_var(ncid, varid , XLONG_V  ))  !"XLONG_V"
!          !call check( nf90_inq_varid(ncid,"XLAT_V" , varid )); call check( nf90_get_var(ncid, varid , XLAT_V   ))  !"XLAT_V"
!
!
!          !lat
!          stat=nf90_def_var(ncid, "lat" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_put_att(ncid, var_id, "units"          , "degrees_north"      )
!          stat=nf90_put_att(ncid, var_id, "long_name"      , "latitude"           )
!          !lon
!          stat=nf90_def_var(ncid, "lon" , NF90_FLOAT, [x_dim_id,y_dim_id], var_id )
!          stat=nf90_put_att(ncid, var_id, "units"          , "degrees_east"       )
!          stat=nf90_put_att(ncid, var_id, "long_name"      , "longitude"          )
!          !var
!          stat=nf90_def_var(ncid, "var" , NF90_FLOAT, [x_dim_id, y_dim_id], var_id )
!          stat=nf90_put_att(ncid, var_id, "units"          , "g/s"                 )
!          stat=nf90_put_att(ncid, var_id, "long_name"      , "var mass flux"       )
!      stat=nf90_enddef(ncid)
!      !Abro NetCDF outFile
!      stat=nf90_open(ncFile, nf90_write, ncid)
!         stat=nf90_inq_varid(ncid, 'var', var_id); stat=nf90_put_var(ncid, var_id, var(:,:)) 
!         stat=nf90_inq_varid(ncid, "lat", var_id); stat=nf90_put_var(ncid, var_id, reshape(lat,[g%nx,g%ny]))!lat(:)) !
!         stat=nf90_inq_varid(ncid, "lon", var_id); stat=nf90_put_var(ncid, var_id, reshape(lon,[g%nx,g%ny]))!lon(:)) !
!      stat=nf90_close(ncid)
!      !Cierro NetCDF outFile
!   end subroutine



!=!================================================================================================
!=!get grid object from diferent files: WRF, CMAQ, GRIDDESC
!= 
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
!=!================================================================================================

end program
