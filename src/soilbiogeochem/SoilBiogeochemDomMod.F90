module SoilBiogeochemDomMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for Dissolved Organic Matter (DOM) carbon production.
  !
  ! !USES:
  use shr_kind_mod                    , only : r8 => shr_kind_r8
  use decompMod                       , only : bounds_type
  !use clm_varcon                      , only : dzsoi_decomp, zisoi
  !TODO add variables clm_varctl
  !use clm_varctl                      , only : use_nitrif_denitrif, use_vertsoilc
  use ColumnType                      , only : col
  use TemperatureType                 , only : temperature_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type                
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: readParams
  public :: SoilBiogeochemDOCprod  ! Calculate DOC production
  !
  ! !PRIVATE DATA:
  type, private :: params_type
     real(r8):: tau_s1_bgc ! 1/turnover time of  SOM 1 from Century (1/7.3) (1/yr)
  end type params_type
  
  type(params_type), private ::  params_inst

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !-----------------------------------------------------------------------

contains
  subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    !
    ! !USES:
    use ncdio_pio    , only: file_desc_t,ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'CNDecompBgcParamsType'
    character(len=100) :: errCode = 'Error reading in CN const file '
    logical            :: readv   ! has variable been read in or not
    real(r8)           :: tempr   ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    !-----------------------------------------------------------------------

    ! Read off of netcdf file
    tString='tau_s1'
    call ncd_io(trim(tString),tempr, 'read', ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%tau_s1_bgc=tempr
  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemDOCprod(bounds, num_soilc, filter_soilc, &
        temperature_inst, totsomc_col, soilbiogeochem_carbonflux_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update the nitrogen leaching rate
    ! as a function of soluble mineral N and total soil water outflow.
    !
    ! !USES:
    use clm_varpar       , only : nlevdecomp, nlevsoi
    use clm_varcon       , only : secspday
    use clm_time_manager , only : get_step_size_real,get_days_per_year
    !
    ! !ARGUMENTS:
    type(bounds_type)                       , intent(in)    :: bounds  
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(temperature_type)                  , intent(in)    :: temperature_inst
    type(soilbiogeochem_carbonflux_type)    , intent(inout) :: soilbiogeochem_carbonflux_inst
    real(r8)                       , intent(in)    :: totsomc_col(bounds%begc:)   ! (gC/m2) total soil organic matter C
    !real(r8)                       , intent(out)   :: somc_doc_col(bounds%begc:)  ! (gC/m2/s) DOC production from soil organic matter

    !SHR_ASSERT_ALL_FL((ubound(totsomc_col)          == (/bounds%endc/)) , sourcefile, __LINE__)
    !SHR_ASSERT_ALL_FL((ubound(somc_dom_col)        == (/bounds%endc/)) , sourcefile, __LINE__) 

    ! !LOCAL VARIABLES:
    integer  :: c,fc                                   ! indices
    real(r8):: kp                                      ! Rate constant for DOC production, specific to
                                                       ! each carbon pool (unit, per day)
    real(r8) :: dt                                     ! radiation time step (seconds)
    real(r8) :: dayspyr                                ! days per year (days)
    real(r8) :: tauz = 24.82748_r8                     ! Empirical factor for a
                                                       ! decrease in C decomposition rates with soil depths
    real(r8), parameter :: Df = 0.75_r8                ! Slope parameter controlling DOC
    real(r8), parameter :: Fs = 1.0_r8                 ! Rate modifier due to soil moisture(kg/m2)
    real(r8) :: Ft                                     ! Rate modifier due to soil temperature(K)

    !-----------------------------------------------------------------------

    associate(                                                                             & 
         t_soisno             =>    temperature_inst%t_soisno_col , & ! Input:  [real(r8) (:,:)  ]  soil temperature (Kelvin) (-nlevsno+1:nlevgrnd) 
         totsomc              =>    totsomc_col, &                    ! Input:  [real(r8) (:)     ]  (gC/m2) total soil organic matter C
         somc_doc             =>    soilbiogeochem_carbonflux_inst%somc_doc_col & ! Output:[real(r8) (:)     ]  (gC/m2/s) DOC production from soil organic matter
         )
      ! set time steps
      dt = get_step_size_real()
      dayspyr = get_days_per_year()
           do fc = 1,num_soilc
              c = filter_soilc(fc)
              kp = 1.e-4_r8    / (secspday * dayspyr * params_inst%tau_s1_bgc)
             !convert Kp unit per day to per seconds
              Ft  = 47.9_r8 / (1 + exp(106._r8 / (t_soisno(c,1) + 18.3_r8)))
              somc_doc(c) = totsomc(c) * (1 - exp ( -kp*Fs*Ft*Df)) * exp(-tauz)
           end do

    end associate

  end subroutine SoilBiogeochemDOCprod

end module SoilBiogeochemDomMod
