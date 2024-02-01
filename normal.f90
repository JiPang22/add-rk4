program ImprovedBoxMullerHistogram
  implicit none

  integer, parameter :: num_samples = 10000000
  real(8), dimension(num_samples) :: random_numbers
  real(8) :: u1, u2, z1, z2
  integer :: i, bins(100), bin_index

  ! 초기화
  bins = 0

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

  ! 히스토그램 계산
  do i = 1, num_samples
    bin_index = min(100, max(1, int(random_numbers(i) * 20.0 + 50.0)))
    bins(bin_index) = bins(bin_index) + 1
  end do
open(1,file='aa')
  ! 히스토그램 출력
  do i = 1, 100
    write(1,*) (i-50)/20.0, bins(i)
  end do

end program ImprovedBoxMullerHistogram

