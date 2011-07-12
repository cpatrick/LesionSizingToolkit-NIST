#
# Insert examples here.
#
option (LSTK_USE_VTK "Build visualization helper tools." OFF)

IF (LSTK_USE_VTK)
  FIND_PACKAGE(VTK REQUIRED)
  INCLUDE(${VTK_USE_FILE})

  SET(ITK_VTK_LSTK_LIBRARIES
   ${ITK-LesionSizingToolkit-Test_LIBRARIES}
   ${ITK-LesionSizingToolkit_LIBRARIES}
   vtkRendering vtkIO vtkWidgets vtkHybrid )

  find_package( Boost REQUIRED )
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

  find_package( LibXml2 )
  if( LIBXML2_FOUND )
    INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
  endif( LIBXML2_FOUND )

 ADD_EXECUTABLE( LesionSegmentation
   LesionSegmentation.cxx
   itkLesionSegmentationCommandLineProgressReporter.cxx )
  TARGET_LINK_LIBRARIES( LesionSegmentation ${ITK_VTK_LSTK_LIBRARIES} )
 
  ADD_EXECUTABLE( LesionSegmentationNIST
   LesionSegmentationNIST.cxx
   itkLesionSegmentationCommandLineProgressReporter.cxx )
  TARGET_LINK_LIBRARIES( LesionSegmentationNIST 
    ${ITK_VTK_LSTK_LIBRARIES} 
    ${LIBXML2_LIBRARIES} )
ENDIF (LSTK_USE_VTK)
