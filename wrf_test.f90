program test

  use netcdf
  use PROJ
  use SCRIP

  implicit none
  character(100), parameter :: proj_latlon = "+proj=latlong +a=6370000.0 +b=6370000.0"
  type(regular_grid_type )  :: g1, g2  !src dst

  !dst: (namelist)
  character(len=256) :: iFile, oFile !src dst
  character(len=20)  :: start_date, end_date 
  integer            :: nx,ny,nz
  real               :: dx,dy,clon,clat
  character(100)     :: proj

  !misc
  real(8)            :: xcent, ycent
  !character(len=256):: method
  integer            :: time,n_hours,soil_levels
  character(10)      :: date_str
  character(8)       :: time_str
  integer :: iostat
  real, dimension(:,:,:,:) :: var1, var2    

  namelist/grid_specs/iFile,oFile, start_date,end_date, nx,ny,nz,soil_levels, dx,dy,clon,clat, proj
  
!---
!(1) Get grid & proj info
   !---------------
   !dst grid specs (from namelist):
   read(*,nml=grid_specs, iostat=iostat)
   if( iostat /= 0 ) stop 'error reading namelist'

   g2%gridName="Dst grid"
   g2%proj4=proj
   g2%nx=nx;  g2%ny=ny;  g2%nz=nz;
   g2%dx=dx;  g2%dy=dy
   xcent=clon; ycent=clat
   call proj_trans(proj_latlon, trim(proj), xcent, ycent)
   g2%xmin=xcent - g2%dx * 0.5 * g2%nx
   g2%ymin=ycent - g2%dy * 0.5 * g2%ny
   g2%xmax=xcent + g2%dx * 0.5 * g2%nx
   g2%ymax=ycent + g2%dy * 0.5 * g2%ny
   !----------------
   !src grids specs:                                              
   date_str=start_date(1:10)
   time_str="00_00_00"!start_date(12:19) 
   iFile=replace(iFile,'<date>', date_str)
   iFile=replace(iFile,'<time>', time_str)
   print*,iFile
   call get_wrf_grid_and_proj(trim(iFile), g1, time, soil_levels)
 
   !debug:
   call print_grid_values(g1)  !debug
   call print_grid_values(g2)  !debug

!---
!(2) Create outFile structure:
   call create_wrf_outfile(oFile, g2)

!---
!(3) Remap variables:
   n_hours=(atoi(date(end_date,"%s"))-atoi(date(start_date,"%s")))/(60*60)
   print*,"# hours:", n_hours
   !allocate(var1(g1%nx,g1%ny))
   !allocate(var2(g2%nx,g2%ny))

   call check(nf90_open(trim(iFile), nf90_read, ncid ))

   do i=1,size(var_list)

       gridDim = grid_dims(i)   ! "2D","3D" (+ time dim)
      gridType = grid_type(i)   ! (M,U,V,W,S)_GRID
       varType =  var_type(i)   ! float, int
       varName =  var_list(i)   !

      g=g1                      !current grid

      select case (trim(gridType))
         case ("U-GRID")
              g%gridName='U'//trim(gridDim)
              g%nx=g%nx+1
              g%xmin=g%xmin-g%dx
              g%xmax=g%xmax+g%dx
         case ("V-GRID")
              g%gridName='V'//trim(gridDim)
              g%ny=g%ny+1
              g%ymin=g%ymin-g%dy
              g%ymax=g%ymax+g%dy
         case ("W-GRID")
              g%gridName='V'//trim(gridDim)
              g%ny=g%ny+1
              g%ymin=g%ymin-g%dy
              g%ymax=g%ymax+g%dy
         case ("S-GRID")
              g%gridName='V'//trim(gridDim)
              g%nz=size(soil_levels)
         case default  !M-GRID
              g%gridName='M'//trim(gridDim)
      end select 

      if (allocated(var1)) deallocate(var1); allocate(var1( g%nx, g%ny,g%nz,Times))  !VER BIEN ESTO
      if (allocated(var2)) deallocate(var2); allocate(var2(g2%nx,g2%ny,g%nz,Times))  !VER BIEN ESTO

      call check( nf90_inq_varid(ncid,trim(varName), varid ))

      do t=1,Times
         select case(type_list(i))
           case ("2D")
             call check(nf90_get_var (ncid, varid, var1(:,:,1,t), start=(/1,1,t/),count=(/nx,ny,1/)) ) 
             call SCRIP_remap_field(var1,var2,g1,g2,method)
           case ("3D")
             do z=1,grid%nz
                call check(nf90_get_var (ncid, varid, var1(:,:,z,t), start=(/1,1,z,t/),count=(/nx,ny,1,1/)) ) 
                call SCRIP_remap_field(var1,var2,g1,g2,method)
             enddo
           case default
             stop 'grid dimension not found.'
         end select
      enddo

      call save_variable_in_netcdf(oFile, varName(i), var2, g2)

      deallocate(tmp)
   enddo

   call check(nf90_close(ncid))

   print*, "Fin de la prueba."

contains

subroutine check(status)
  integer, intent(in) :: status
  if (status /= nf90_noerr) then
    write(*,*) nf90_strerror(status)
    stop 'netcdf error'
  end if
end subroutine check
!Utils: ----------------------------------------
integer function atoi(str)            !string -> int
  implicit none
  character(len=*), intent(in) :: str
  read(str,*) atoi
end function
character(len=20) function itoa(i)    !int -> string
   implicit none
   integer, intent(in) :: i
   write(itoa, '(i0)') i
   itoa = adjustl(itoa)
end function
character(len=16) function rtoa(r)    !real -> string
   implicit none
   real, intent(in) :: r
   write(rtoa, '(F16.3)') r
   rtoa = adjustl(rtoa)
end function
pure function replace(string, s1, s2) result(str)
    !replace substring "s1" with "s2" on string
    implicit none
    character(*), intent(in)       :: string
    character(*), intent(in)       :: s1,s2
    character(len(string)+len(s2)) :: str          !not very elegant
    integer :: i,j,n,n1,n2!,dif
    str=string;n =len(str); n1=len(s1); n2=len(s2)
    if ( n1 == n2 ) then
       do i=1,n-n1
          if ( str(i:i+n1) == s1 ) str(i:i+n1) = s2
       end do
    else
        i=1
        do while ( i < len(trim(str)))
           if ( str(i:i+n1-1) == s1 ) then
                str(i+n2:n)=str(i+n1:n)            !make space on str for replacement
                str(i:i+n2-1) = s2                 !replace! 
           endif
          i=i+1
        enddo
    endif
end function

character(len=20) function date(date_str, fmt_str) result(output) !Interfaz a "date"
  implicit none
  integer :: iostat
  character(*), intent(in) :: date_str, fmt_str
  character(256)           :: command
  command="date -d '"//trim(date_str)//"' '+"//trim(fmt_str)//"'  > tmp.date"
  call system( trim(command) )
  open(9, file='tmp.date', status='old',action='read'); read(9, '(A)', iostat=iostat) output;  close(9)
  call system('rm tmp.date')
end function
!-----------------------------------------------

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
    g%gridName='Src grid'

    call check(nf90_open(trim(wrf_file), nf90_write, ncid ))
       !proj params:                    
       call check(nf90_get_att(ncid, NF90_GLOBAL, "MAP_PROJ",  proj_id  )) !1: lcc  2: stere 3: merc  6: latlon
       call check(nf90_get_att(ncid, NF90_GLOBAL, "CEN_LAT" ,  CEN_LAT  ))
       call check(nf90_get_att(ncid, NF90_GLOBAL, "CEN_LON" ,  CEN_LON  ))
       call check(nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT1",  TRUELAT1 ))
       call check(nf90_get_att(ncid, NF90_GLOBAL, "TRUELAT2",  TRUELAT2 ))
       call check(nf90_get_att(ncid, NF90_GLOBAL, "STAND_LON", STAND_LON))
       !grid params:    
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
      write(g%proj4,'("+proj=lcc +a=6370000.0 +b=6370000.0 +lat_1=",F8.4," +lat_2=",F8.4," +lat_0=",F8.4," +lon_0=",F8.4)') truelat1,truelat2,cen_lat,cen_lon!stand_lon!
     case (2) !stere
      write(g%proj4,'("+proj=stere +a=6370000.0 +b=6370000.0 +lat_0=",F9.4," +lon_0=",F9.4," +lat_ts=",F9.4)') truelat1,truelat2,cen_lat,stand_lon
     case (3) !merc
      write(g%proj4,'("+proj=merc +a=6370000.0 +b=6370000.0 +lon_0=",F9.4," +lat_ts=",F9.4)') truelat1,truelat2,stand_lon,cen_lat
     case (6)
      write(g%proj4,'("+proj=latlong +a=6370000.0 +b=6370000.0")')
    end select

    !set Grid parameters:
    xcent=CEN_LON; ycent=CEN_LAT
    call proj_trans(proj_latlon, g%proj4, xcent, ycent )
    !M-GRID
    g%xmin=xcent - g%dx *0.5* g%nx
    g%xmax=xcent + g%dx *0.5* g%nx
    g%ymin=ycent - g%dy *0.5* g%ny
    g%ymax=ycent + g%dy *0.5* g%ny
end subroutine


   subroutine create_wrf_outfile(ncFile,g2)
      implicit none
      character(len=*),intent(in) :: ncfile
      type(regular_grid_type )    :: g1,g2  !src dst
      !Coordinates:
      real(8), allocatable,dimension(:,:) :: XLAT,XLAT_U,XLAT_V
      real(8), allocatable,dimension(:,:) :: XLONG,XLONG_U,XLONG_V
      real(8), allocatable,dimension(:)   :: Times,ZNU,ZNW,ZS     
      integer :: i,j,stat,idx
      character(100) :: proj4_latlon="+proj=latlong +a=6370000.0 +b=6370000.0"

      integer :: ncid,varid
      integer :: dim_dateStr,dim_Time,dim_x,dim_y,dim_z ,dim_x_stag ,dim_y_stag ,dim_z_stag ,dim_soil 

      integer :: soilLayers=4
      allocate(  Times( 24)              ) !daily files
      allocate(   XLAT( g2%nx  , g2%ny  ))
      allocate( XLAT_U( g2%nx+1, g2%ny  ))
      allocate( XLAT_V( g2%nx  , g2%ny+1))
      allocate(  XLONG( g2%nx  , g2%ny  ))
      allocate(XLONG_U( g2%nx+1, g2%ny  ))
      allocate(XLONG_V( g2%nx  , g2%ny+1))
      allocate(    ZNU( g2%nz     )      )
      allocate(    ZNW( g2%nz+1   )      )
      allocate(     ZS( soilLayers))

      ! !1D coords:
      ! Times=t  !con start/end_date
      !   ZNU=   !aprender como calcular en base a ptop y n_vert!
      !   ZNW=
      !    ZS=   !supongo que depende del Land Surface Model
      !2d coords:
      do j=1,g2%ny
        do i=1,g2%nx
            XLONG(i,j)=g2%xmin+0.5*g2%dx+g2%dx*i
             XLAT(i,j)=g2%ymin+0.5*g2%dy+g2%dy*j
        enddo
      enddo                                                                   
      
      XLONG_U(1:g2%nx,1:g2%ny)=XLONG - g2%dx; XLONG_U(g2%nx+1,:)=XLONG(g2%nx,:)+g2%dx  
       XLAT_U(1:g2%nx,1:g2%ny)=XLAT         ;  XLAT_U(g2%nx+1,:)= XLAT(g2%nx,:)       
      XLONG_V(1:g2%nx,1:g2%ny)=XLONG        ; XLONG_V(:,g2%ny+1)=XLONG(:,g2%ny)       
       XLAT_V(1:g2%nx,1:g2%ny)=XLAT  - g2%dy;  XLAT_V(:,g2%ny+1)= XLAT(:,g2%ny)+g2%dy  

      call proj_trans(g2%proj4, proj4_latlon, Xlong  , Xlat  )
      call proj_trans(g2%proj4, proj4_latlon, Xlong_U, Xlat_U)  !Ineficiente
      call proj_trans(g2%proj4, proj4_latlon, Xlong_V, Xlat_V)  !Ineficiente
      !print*,"minval lat,lon:",minval(xlat),minval(xlong)
      !print*,"maxval lat,lon:",maxval(xlat),maxval(xlong)

      print*,"write NETCDF"
      stat=nf90_create(ncFile, NF90_CLOBBER, ncid)
          !dims
          stat=nf90_def_dim(ncid, "DateStrLen"      , 19        , dim_dateStr)
          stat=nf90_def_dim(ncid, "Time"            , 24        , dim_Time   ) !daily files
          stat=nf90_def_dim(ncid, "west_east"       , g2%nx     , dim_x      )
          stat=nf90_def_dim(ncid, "south_north"     , g2%ny     , dim_y      )
          stat=nf90_def_dim(ncid, "bottom_top"      , g2%nz     , dim_z      )
          stat=nf90_def_dim(ncid, "west_east_stag"  , g2%nx+1   , dim_x_stag )
          stat=nf90_def_dim(ncid, "south_north_stag", g2%ny+1   , dim_y_stag )
          stat=nf90_def_dim(ncid, "bottom_top_stag" , g2%nz+1   , dim_z_stag )
          stat=nf90_def_dim(ncid, "soil_layers_stag", soilLayers, dim_soil   )

          !def coord vars:
          stat=nf90_def_var(ncid, "Times"  ,  NF90_FLOAT, [dim_dateStr,dim_Time], varId )
          stat=nf90_def_var(ncid, "ZS"     ,  NF90_FLOAT, [dim_soil  ]          , varId ) 
          stat=nf90_def_var(ncid, "ZNU"    ,  NF90_FLOAT, [dim_z     ]          , varId )
          stat=nf90_def_var(ncid, "ZNW"    ,  NF90_FLOAT, [dim_z_stag]          , varId )
          stat=nf90_def_var(ncid, "XLONG"  ,  NF90_FLOAT, [dim_x, dim_y]        , varId )
          stat=nf90_def_var(ncid, "XLAT"   ,  NF90_FLOAT, [dim_x, dim_y]        , varId )
          stat=nf90_def_var(ncid, "XLONG_U",  NF90_FLOAT, [dim_x_stag, dim_y]   , varId )
          stat=nf90_def_var(ncid, "XLAT_U" ,  NF90_FLOAT, [dim_x_stag, dim_y]   , varId )
          stat=nf90_def_var(ncid, "XLONG_V",  NF90_FLOAT, [dim_x, dim_y_stag]   , varId )
          stat=nf90_def_var(ncid, "XLAT_V" ,  NF90_FLOAT, [dim_x, dim_y_stag]   , varId )
          !stat=nf90_put_att(ncid, var_id, "units"          , "g/s"                 )
          !stat=nf90_put_att(ncid, var_id, "long_name"      , "var mass flux"       )
      stat=nf90_enddef(ncid)

      !Add variables to netcdf:
      stat=nf90_open(ncFile, nf90_write, ncid)
         !stat=nf90_inq_varid(ncid, "Times"  , varId); stat=nf90_put_var(ncid, var_id,        ) 
         !stat=nf90_inq_varid(ncid, "ZS"     , varId); stat=nf90_put_var(ncid, var_id,        ) 
         !stat=nf90_inq_varid(ncid, "ZNU"    , varId); stat=nf90_put_var(ncid, var_id,        ) 
         !stat=nf90_inq_varid(ncid, "ZNW"    , varId); stat=nf90_put_var(ncid, var_id,        ) 
         stat=nf90_inq_varid(ncid, "XLONG"  , varId); stat=nf90_put_var(ncid, varId, XLONG   ) 
         stat=nf90_inq_varid(ncid, "XLAT"   , varId); stat=nf90_put_var(ncid, varId, XLAT    ) 
         stat=nf90_inq_varid(ncid, "XLONG_U", varId); stat=nf90_put_var(ncid, varId, XLONG_U ) 
         stat=nf90_inq_varid(ncid, "XLAT_U" , varId); stat=nf90_put_var(ncid, varId, XLAT_U  ) 
         stat=nf90_inq_varid(ncid, "XLONG_V", varId); stat=nf90_put_var(ncid, varId, XLONG_V ) 
         stat=nf90_inq_varid(ncid, "XLAT_V" , varId); stat=nf90_put_var(ncid, varId, XLAT_V  ) 
      stat=nf90_close(ncid)
      !Cierro NetCDF outFile
   end subroutine

!degug: ----------------------------------------
subroutine print_grid_values(g)
  implicit none
  type(regular_grid_type) , intent(in) :: g
  real (8) :: minlat,maxlat,minlon,maxlon

  minlat=g%ymin;maxlat=g%ymax
  minlon=g%xmin;maxlon=g%xmax
  call proj_trans(g%proj4, proj_latlon,minlon,minlat ) 
  call proj_trans(g%proj4, proj_latlon,maxlon, maxlat) 
   !debug:
   print '("---",/,A10)', g%gridName
   print '("proj:", A110)', g%proj4
   print '("nx:",I4," ny:",I4," nz:", I4," soil_layers:", I4)',g%nx, g%ny, g%nz, 4
   print '("dx, dy:",2(F9.3))' ,g%dx, g%dy
   print '("  x/y   min:",2(F12.3),/,"  x/y   max:",2(F12.3))', g%xmin,g%ymin,g%xmax, g%ymax
   print '("lon/lat-min:",2(F12.3),/,"lon/lat-max:",2(F12.3))', minlon,minlat,maxlon, maxlat
end subroutine

end program
