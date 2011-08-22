PROJECT( HOUGHCIRCLESEGMENTATION )

# CMake >= 2.6 MUST be installed
CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )

IF( COMMAND CMAKE_POLICY )
  CMAKE_POLICY( VERSION 2.4 )
  CMAKE_POLICY( SET CMP0005 OLD )
  CMAKE_POLICY( SET CMP0003 NEW )
ENDIF(COMMAND CMAKE_POLICY)

#---------------------------------------------------------------------------
# Find ITK.
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR
    "Cannot build without ITK.  Please set ITK_DIR.")
ENDIF(ITK_FOUND)

INCLUDE( CTest )

ADD_EXECUTABLE(hough2D HoughTransformRadialVotingImageFilter2D.cxx)
TARGET_LINK_LIBRARIES(hough2D ITKBasicFilters ITKAlgorithms ITKStatistics ITKNumerics ITKCommon ITKIO)

ADD_EXECUTABLE(hough3D HoughTransformRadialVotingImageFilter3D.cxx)
TARGET_LINK_LIBRARIES(hough3D ITKBasicFilters ITKAlgorithms ITKStatistics ITKNumerics ITKCommon ITKIO)

ADD_TEST( hough2DTest
  ${EXECUBTALE_OUTPUT_PATH}/hough2D
  ${HOUGHCIRCLESEGMENTATION_SOURCE_DIR}/input.tif
  output.png
  accumulator.mha
  50 6 7.5 1 1 1 0.15 100 0.5 1 1 1
)