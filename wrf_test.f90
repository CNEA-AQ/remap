program test

  use netcdf
  use proj  
  !use remap
  use SCRIP

  implicit none
  type(regular_grid_type ) :: g1,g2
  character(len=256) :: srcFile, dstFile, outFile
  character(len=256) :: method
  real, allocatable :: var1(:,:),var2(:,:)
  
  type wrf_file
    character(*) :: fileName
    type(grid_type) :: m_grd, u_grd, v_grd, w_grd, s_grd !horizontal and vertical grids
  end type
  type (wrf_file) :: wrf

   call system('echo "(!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!)"');  
   call system('echo "(!) BORRO TODOS LOS NETCDF EN DIRECTORIO    "');  
   call system('echo "                                            "');  
   call system('rm *.nc');  

   srcFile='/media/usuario/cnea_ram/WRFE/2023112400/01/wrfout_d01_2023-11-24_00_00_00' !'wrfout_d01_20XX-XX-XX_XX:00:00.nc'
   dstFile='./template/wrfout_base.nc'
   dstFile='./wrfout_remapeado.nc'

   !***************************
   !src grids specs:
   call get_wrf_wrid_and_proj(srcFile, g1)
   allocate(var1(g1%nx,g1%ny))
   !***************************
   !dst grid specs (silam grid):
   call get_wrf_grid_and_proj(dstFile, g2)
   allocate(var2(g2%nx,g2%ny))

   !***************************
   !Remap variables:
   do i=1,size(var_list)

      gridType = grid_type(i)   ! (M,U,V,W,S)_GRID
       varType =  var_type(i)   ! float, int

      allocate(tmp(Times,nx,ny,nz))

      do t=1,wrf%Times
        do z=1,grid%nz

           select case(type_list(i))
             case ("")
                  call SCRIP_remap_field(var1,var2,g1,g2,method)
             case ("")
             case ("")
             case ("")
        enddo
      enddo

      !Save variable in outFile:

      deallocate(tmp)
   enddo

   !!!test Bilinear   
   !method='bilinear'
   !ofile='ou_'//trim(method)//'.nc'
   !call SCRIP_remap_field(var1,var2,g1,g2,method)
   !print*, "Bilinear Remaping succesfull!"
   !call saveArrayOnNetCDF(oFile,g2,var2)

   call saveArrayOnNetCDF(outFile,g2,var2)


   print*, "Fin de la prueba."

contains

   subroutine get_wrf_grid_and_proj(g,var)
       implicit none
       type(regular_grid_type ) :: g
       real                :: var(:,:)
       integer :: i,j,stat!,nx,ny
       !proj
       real    :: cen_lat, cen_lon, truelat1,truelat2,stand_lon
       real    :: a=6370000.0,b=6370000.0      !axis and semi-axis radius of wrf-spheroid
       character(len=*) :: proj_str
       !grid
       integer :: south_north,west_east,bottom_top,time
       integer :: dx,dy 

       call check(nf90_open(trim(inp_file), nf90_write, ncid ))
            !Projection parameters
            call check(nf90_get_att(ncid, NF90_GLOBAL, "MAP_PROJ",  proj_id  )) !1: lcc  2: stere 3: merc  6: latlon
            call check(nf90_get_att(ncid, NF90_GLOBAL, "CEN_LAT" ,  CEN_LAT  ))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "CEN_LON" ,  CEN_LON  ))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1",  TRUELAT1 ))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2",  TRUELAT2 ))
            call check(nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", STAND_LON))

            !Grid dimensions:
            !DateStrLen          = 19;
            call check(nf90_inq_dimid(ncid,"Time"            ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%Time            ) )
            call check(nf90_inq_dimid(ncid,"south_north"     ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%south_north     ) )
            call check(nf90_inq_dimid(ncid,"west_east"       ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%west_east       ) )
            call check(nf90_inq_dimid(ncid,"bottom_top"      ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%bottom_top      ) )
            call check(nf90_inq_dimid(ncid,"west_east_stag"  ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%west_east_stag  ) )
            call check(nf90_inq_dimid(ncid,"south_north_stag",dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%south_north_stag) )
            call check(nf90_inq_dimid(ncid,"bottom_top_stag" ,dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%bottom_top_stag ) )
            call check(nf90_inq_dimid(ncid,"soil_layers_stag",dimid)); call check(nf90_inquire_dimension(ncid,dimid,len=wf%soil_layers_stag) )
       
            call check(nf90_get_att(ncid, NF90_GLOBAL, "DX", grid%dx) )
            call check(nf90_get_att(ncid, NF90_GLOBAL, "DY", grid%dy) )                                                                                                                  
       call check(nf90_close(ncid))
     
       !set Projection       
       select case (proj_id)
         case (1) !lcc
            ! "+proj=lcc     +a=6370000.0 +b=6370000.0 +lat_1=33 +lat_2=45 +lat_0=40 +lon_0=-97", TRUELAT1,TRUELAT2,CEN_LAT,CEN_LON
            write(proj_str,'("+proj=lcc +a=",F9.4," +b=",F9.4," +lat_1=",F9.4," +lat_2=",F9.4," +lat_0=",F9.4," +lon_0=",F9.4)') a,b,truelat1,truelat2,cen_lat,stand_lon
         case (2) !stere
            ! "+proj=stere   +a=6370000.0 +b=6370000.0 +lat_ts=33 +lat_0=90 +lon_0=-97 +k_0=1.0"
            write(proj_str,'("+proj=stere +a=",F9.4," +b=",F9.4," +lat_0=",F9.4," +lon_0=",F9.4," +lat_ts=",F9.4)') a,b,truelat1,truelat2,cen_lat,stand_lon
         case (3) !merc
            ! "+proj=merc    +a=6370000.0 +b=6370000.0 +lat_ts=33 +lon_0=0"
            write(proj_str,'("+proj=merc +a=",F9.4," +b=",F9.4," +lon_0=",F9.4," +lat_ts=",F9.4)') a,b,truelat1,truelat2,stand_lon,cen_lat
         case (6)
           ! "+proj=latlong +a=6370000.0 +b=6370000.0"
           write(proj_str,'("+proj=latlong +a=",F9.1," b=",F9.1"') a,b
       end select

       wrf%m_grid%proj=proj_str
       wrf%u_grid%proj=proj_str
       wrf%v_grid%proj=proj_str
       wrf%w_grid%proj=proj_str
       wrf%w_grid%proj=proj_str

       !set Grid parameters:
       xcent=CEN_LON; ycent=CEN_LAT
       call proj_trans(proj_latlon, proj_str, xcent,ycent)

       !M-GRID
       grid%xmin=xcent - grid%dx * grid%nx
       grid%ymin=ycent - grid%ny * grid%ny
       grid%xmax=xcent + grid%dx * grid%nx
       grid%ymax=ycent + grid%ny * grid%ny
       
       !U-GRID

       !V-GRID

       !W-GRID

       !M-GRID

       !!Coordinates:
       !allocate( Times(Time)) ;
       !allocate(    XLAT(south_north, west_east     )) ;!Time,
       !allocate(  XLAT_U(south_north, west_east_stag)) ;!Time,
       !allocate(  XLAT_V(south_north_stag, west_east)) ;!Time,
       !allocate(   XLONG(south_north, west_east     )) ;!Time,
       !allocate( XLONG_U(south_north, west_east_stag)) ;!Time,
       !allocate( XLONG_V(south_north_stag, west_east)) ;!Time,
       !allocate(     ZNU(bottom_top                 )) ;!Time,
       !allocate(     ZNW(bottom_top_stag            )) ;!Time,
       !allocate(      ZS(soil_layers_stag           )) ;!Time,
                                                                                                                            
       !call check( nf90_inq_varid(ncid,"Times"  , varid )); call check( nf90_get_var(ncid, varid , Times    ))  !"Times"
       !call check( nf90_inq_varid(ncid,"ZS"     , varid )); call check( nf90_get_var(ncid, varid , ZS       ))  !"ZS"
       !call check( nf90_inq_varid(ncid,"ZNU"    , varid )); call check( nf90_get_var(ncid, varid , ZNU      ))  !"ZNU"
       !call check( nf90_inq_varid(ncid,"ZNW"    , varid )); call check( nf90_get_var(ncid, varid , ZNW      ))  !"ZNW"
       !call check( nf90_inq_varid(ncid,"XLONG"  , varid )); call check( nf90_get_var(ncid, varid , XLONG    ))  !"XLONG"
       !call check( nf90_inq_varid(ncid,"XLAT"   , varid )); call check( nf90_get_var(ncid, varid , XLAT     ))  !"XLAT"
       !call check( nf90_inq_varid(ncid,"XLONG_U", varid )); call check( nf90_get_var(ncid, varid , XLONG_U  ))  !"XLONG_U"
       !call check( nf90_inq_varid(ncid,"XLAT_U" , varid )); call check( nf90_get_var(ncid, varid , XLAT_U   ))  !"XLAT_U"
       !call check( nf90_inq_varid(ncid,"XLONG_V", varid )); call check( nf90_get_var(ncid, varid , XLONG_V  ))  !"XLONG_V"
       !call check( nf90_inq_varid(ncid,"XLAT_V" , varid )); call check( nf90_get_var(ncid, varid , XLAT_V   ))  !"XLAT_V"

   end subroutine


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
