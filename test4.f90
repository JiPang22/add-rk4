program BoxMullerGaussianTest
  use, intrinsic :: iso_c_binding
  use :: gsl_types
  implicit none

  integer, parameter :: num_samples = 10000
  real(8), dimension(num_samples) :: random_numbers
  real(8) :: u1, u2, z1, z2
  integer :: i, status
  type(c_ptr) :: ws
  real(8) :: wstat, pvalue

  ! GSL 샤피로-윌크 검정을 위한 변수
  integer, parameter :: n = num_samples
  real(8), dimension(n) :: data

  call random_seed()  ! 시드 설정
  do i = 1, num_samples / 2
    ! 두 개의 독립적인 균등분포 난수 생성
    call random_number(u1)
    call random_number(u2)

    ! Box-Muller 변환 수행
    z1 = sqrt(-2.0d0 * log(u1)) * cos(2.0d0 * 3.141592653589793d0 * u2)
    z2 = sqrt(-2.0d0 * log(u1)) * sin(2.0d0 * 3.141592653589793d0 * u2)

    ! 생성된 정규 분포 난수 저장
    random_numbers(2 * i - 1) = z1
    random_numbers(2 * i) = z2
  end do

  ! GSL 샤피로-윌크 검정을 위해 데이터 복사
  data = random_numbers

  ! GSL 초기화
  call gsl_set_error_handler(0_cptr)
  ws = gsl_stats_wstat_alloc(n)
  status = gsl_stats_wshapiro_test(data, 1_csize_t, ws)

  ! 샤피로-윌크 검정 통계량 및 p-value 획득
  wstat = gsl_stats_wstat_wkurt(ws)
  pvalue = gsl_stats_wstat_p(ws)

  ! 결과 출력
!  print *, "Shapiro-Wilk test statistic:", wstat
  print *, "P-value:", pvalue

  ! 메모리 해제
  call gsl_stats_wstat_free(ws)

end program BoxMullerGaussianTest

