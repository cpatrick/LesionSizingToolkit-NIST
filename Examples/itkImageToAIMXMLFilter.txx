/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkImageToAIMXMLFilter.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkImageToAIMXMLFilter_txx
#define _itkImageToAIMXMLFilter_txx

#include <time.h>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include <libxml/parser.h>
#include <libxml/tree.h>

#include "itkImageToAIMXMLFilter.h"

namespace itk
{

/**
 * Constructor
 */
template <class TInputImage, class TReferenceImage>
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::ImageToAIMXMLFilter()
  :  m_InputImage(NULL),
     m_ReferenceImage(NULL),
     m_ContourThreshold(0),
     m_Output(""),
     m_CurrentUID(""),
     m_PatientName(""),
     m_PatientId(""),
     m_PatientSex(""),
     m_StudyInstanceUID(""),
     m_SeriesInstanceUID(""),
     m_SOPClassUIDs(),
     m_SOPInstanceUIDs(),
     m_Contours()
{
}

/**
 * Destructor
 */
template <class TInputImage, class TReferenceImage>
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::~ImageToAIMXMLFilter()
{
}

/**
 * Set an itk::Image as input 
 */
template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::SetInput( const InputImageType * inputImage )
{
  m_InputImage = inputImage;
}

/**
 * Set an itk::Image as input 
 */
template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::SetReference( const ReferenceImageType * referenceImage )
{
  m_ReferenceImage = referenceImage;
}

template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::SetSOPClassUIDs( const UIDContainerType& uids )
{
  m_SOPClassUIDs = uids;
}

template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::SetSOPInstanceUIDs( const UIDContainerType& uids )
{
  m_SOPInstanceUIDs = uids;
}

//----------------------------------------------------------------------------
template <class TInputImage, class TReferenceImage>
void 
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::AddContourToVector( typename ContourFilterType::VertexListConstPointer verts,
                      unsigned int sliceNum,
                      ContourPointVectorType& vect,
                      ReferenceIndexValueType& index ) const
{
  // We keep referenceIndex around to populate the key, index, when the loop
  // completes
  typename ReferenceImageType::IndexType referenceIndex;

  // Iterate through the vertices and add them to the output index vector,
  // vect. The transformations move the index from the input images space into
  // world space and the into the reference image's index space.
  for( unsigned int vertexNum = 0; vertexNum < verts->Size(); ++vertexNum )
    {
    typename ContourFilterType::VertexType vertex;
    typename InputImageType::IndexType inputIndex;
    typename InputImageType::PointType point;
    vertex = verts->ElementAt(vertexNum);
    inputIndex[0] = vertex[0];
    inputIndex[1] = vertex[1];
    inputIndex[2] = sliceNum;
    m_InputImage->TransformIndexToPhysicalPoint( inputIndex, point );
    m_ReferenceImage->TransformPhysicalPointToIndex( point, referenceIndex);
    vect.push_back(referenceIndex);
    }
  index = referenceIndex[2];
}

//----------------------------------------------------------------------------
template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::PopulateContourContainer()
{
  // Setup the slice extractor and regions
  typename ExtractFilterType::Pointer extractor = ExtractFilterType::New();
  typename InputImageType::RegionType inputRegion =
    m_InputImage->GetLargestPossibleRegion();
  typename InputImageType::SizeType inputSize = inputRegion.GetSize();
  typename InputImageType::IndexType inputIndex = inputRegion.GetIndex();
  typename InputImageType::SizeType desiredSize = inputSize;
  desiredSize[2] = 0;
  extractor->SetInput( m_InputImage );

  // Iterate through the slices and determine contours using our contour
  // extractor
  for( unsigned int sliceNum = inputIndex[2]; 
       sliceNum < inputSize[2]; 
       ++sliceNum )
    {
    typename InputImageType::RegionType desiredRegion;
    typename InputImageType::IndexType desiredIndex = inputIndex;
    typename ContourFilterType::Pointer contourer = ContourFilterType::New();
    inputIndex[2] = sliceNum;
    desiredRegion.SetSize(  desiredSize  );
    desiredRegion.SetIndex( desiredIndex );

    extractor->SetExtractionRegion( desiredRegion );
    extractor->SetDirectionCollapseToIdentity();

    contourer->SetInput( extractor->GetOutput() );
    contourer->SetContourValue( m_ContourThreshold );
    contourer->Update();

    // The contourer may generate multiple contours per slice, so we iterate
    // through them and add the vectors of them to the multimap of index
    // vectors
    for( unsigned int contourNum = 0;
         contourNum < contourer->GetNumberOfOutputs();
         ++contourNum )
      {
      typename ContourFilterType::VertexListConstPointer verts;
      verts = contourer->GetOutput(contourNum)->GetVertexList();
      ContourPointVectorType contourVector;
      IndexValueType referenceSlice;
      this->AddContourToVector( verts, sliceNum, contourVector, referenceSlice );
      typename ContourContainerType::value_type pairForInsertion(referenceSlice,
                                                                 contourVector);
      m_Contours.insert(pairForInsertion);
      }
    }
}

//----------------------------------------------------------------------------
template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::GenerateXMLFromContours()
{

  // Constants that seem to be inherent to the format, yet non-variable
  const char* XmlVersion = "1.0";
  const char* Encoding = "UTF-8";
  const char* XMLNS =
    "gme://caCORE.caCORE/3.2/edu.northwestern.radiology.AIM"; 
  const char* XMLNSXSI =
    "http://www.w3.org/2001/XMLSchema-instance";
  const char* CodingSchemeVersion = "v0_rv1";
  const char* CodingSchemeDesignator = "AVT2";
  const char* CodeMeaning = "Tumor volume contours";
  const char* CodeValue = "AVT006";
  const char* CompanyName = "SIEMENS";
  const char* AimVersion = "TCGA";
  const char* XSISchemaLocation = 
    "gme://caCORE.caCORE/3.2/edu.northwestern.radiology.AIM "
    "AIM_TCGA09302009_XML.xsd";
  const char* DefaultNumber = "0";
  const char* DefaultName = "LSTK";

  // Creating an XML document using libxml2
  xmlDocPtr doc = NULL;
  xmlNodePtr rootNode;
  xmlNodePtr curNode;
  xmlChar* xmlBuffer;
  int xmlBufferSize;
  unsigned int count;
  ReferenceIndexValueType lastValue;

  // Create new document with requisite properties
  doc = xmlNewDoc(BAD_CAST XmlVersion);
  doc->standalone = 1;
  doc->encoding = xmlCharStrdup( Encoding );

  // Create the big root node required for image annotations and add it
  // to the document as the root.
  rootNode = xmlNewNode(NULL, BAD_CAST "ImageAnnotation");
  xmlNewProp( rootNode, BAD_CAST "xmlns", BAD_CAST XMLNS );
  xmlNewProp( rootNode, BAD_CAST "xmlns:xsi", BAD_CAST XMLNSXSI );
  xmlNewProp( rootNode, 
              BAD_CAST "codingSchemeVersion",
              BAD_CAST CodingSchemeVersion );
  xmlNewProp( rootNode, 
              BAD_CAST "codingSchemeDesignator",
              BAD_CAST CodingSchemeDesignator );
  xmlNewProp( rootNode, 
              BAD_CAST "codeMeaning",
              BAD_CAST CodeMeaning );
  xmlNewProp( rootNode, 
              BAD_CAST "codeValue",
              BAD_CAST CodeValue );
  xmlNewProp( rootNode, 
              BAD_CAST "uniqueIdentifier",
              BAD_CAST m_CurrentUID.c_str() );
  xmlNewProp( rootNode, 
              BAD_CAST "name",
              BAD_CAST CompanyName );
  xmlNewProp( rootNode, 
              BAD_CAST "dateTime",
              BAD_CAST this->GetCurrentTimeInAIMFormat().c_str() );
  xmlNewProp( rootNode, 
              BAD_CAST "aimVersion",
              BAD_CAST AimVersion );
  xmlNewProp( rootNode, 
              BAD_CAST "xsi:schemaLocation",
              BAD_CAST XSISchemaLocation );
  xmlDocSetRootElement(doc, rootNode);

  // Add the user->User tags
  curNode = xmlNewChild(rootNode, NULL, BAD_CAST "user", NULL);
  curNode = xmlNewChild(curNode, NULL, BAD_CAST "User", NULL);
  xmlNewProp( curNode, 
              BAD_CAST "numberWithinRoleOfClinicalTrial",
              BAD_CAST DefaultNumber );
  xmlNewProp( curNode, 
              BAD_CAST "roleInTrial",
              BAD_CAST "Referring" );
  xmlNewProp( curNode, 
              BAD_CAST "loginName",
              BAD_CAST DefaultName );
  xmlNewProp( curNode, 
              BAD_CAST "name",
              BAD_CAST DefaultName );
  xmlNewProp( curNode, 
              BAD_CAST "id",
              BAD_CAST DefaultNumber );

  // Add image reference collection tags
  curNode = xmlNewChild(rootNode, NULL, 
                        BAD_CAST "imageReferenceCollection", NULL);
  curNode = xmlNewChild(curNode, NULL, BAD_CAST "ImageReference", NULL);
  xmlNewProp( curNode,
              BAD_CAST "xsi:type",
              BAD_CAST "DICOMImageReference" );
  xmlNewProp( curNode,
              BAD_CAST "id",
              BAD_CAST DefaultNumber );
  curNode = xmlNewChild(curNode, NULL, BAD_CAST "study", NULL);
  curNode = xmlNewChild(curNode, NULL, BAD_CAST "Study", NULL);
  xmlNewProp( curNode,
              BAD_CAST "instanceUID",
              BAD_CAST m_StudyInstanceUID.c_str() );
  xmlNewProp( curNode,
              BAD_CAST "id",
              BAD_CAST DefaultNumber );
  curNode = xmlNewChild(curNode, NULL, BAD_CAST "series", NULL);
  curNode = xmlNewChild(curNode, NULL, BAD_CAST "Series", NULL);
  xmlNewProp( curNode,
              BAD_CAST "instanceUID",
              BAD_CAST m_SeriesInstanceUID.c_str() );
  xmlNewProp( curNode,
              BAD_CAST "id",
              BAD_CAST DefaultNumber );
  curNode = xmlNewChild( curNode, NULL, BAD_CAST "imageCollection", NULL);

  // Iterate through the multimap to pull out the reference image UIDs
  count = 0;
  lastValue = -1; // worried about this
  for( typename ContourContainerType::const_iterator itr = m_Contours.begin();
       itr != m_Contours.end();
       ++itr )
    {
    if( itr->first != lastValue )
      {
      std::string imageId = boost::lexical_cast<std::string>(count);
      xmlNodePtr imageNode = xmlNewChild( curNode, NULL, 
                                          BAD_CAST "Image", NULL );
      xmlNewProp( imageNode,
                  BAD_CAST "sopInstanceUID",
                  BAD_CAST m_SOPInstanceUIDs.at(itr->first).c_str() );
      xmlNewProp( imageNode,
                  BAD_CAST "sopClassUID",
                  BAD_CAST m_SOPClassUIDs.at(itr->first).c_str() );
      xmlNewProp( imageNode,
                  BAD_CAST "id",
                  BAD_CAST imageId.c_str() );
      ++count;
      lastValue = itr->first;
      }
    }

  // Setup the patient->Patient tags
  curNode = xmlNewChild(rootNode, NULL, BAD_CAST "patient", NULL);
  curNode = xmlNewChild(curNode, NULL, BAD_CAST "Patient", NULL);
  xmlNewProp( curNode,
              BAD_CAST "sex",
              BAD_CAST m_PatientSex.c_str() );
  xmlNewProp( curNode,
              BAD_CAST "patientID",
              BAD_CAST m_PatientId.c_str() );
  xmlNewProp( curNode,
              BAD_CAST "name",
              BAD_CAST m_PatientName.c_str() );
  xmlNewProp( curNode,
              BAD_CAST "id",
              BAD_CAST DefaultNumber );

  // Get Started with the geometric shape collection
  curNode = xmlNewChild( rootNode, NULL,
                         BAD_CAST "geometricShapeCollection", NULL );
  curNode = xmlNewChild( curNode, NULL,
                         BAD_CAST "GeometricShape", NULL );
  xmlNewProp( curNode,
              BAD_CAST "xsi:type",
              BAD_CAST "MultiPoint" );
  xmlNewProp( curNode,
              BAD_CAST "shapeIdentifier",
              BAD_CAST "1" );
  xmlNewProp( curNode,
              BAD_CAST "includeFlag",
              BAD_CAST "true" );
  xmlNewProp( curNode,
              BAD_CAST "lineStyle",
              BAD_CAST "SOLID" );
  xmlNewProp( curNode,
              BAD_CAST "lineOpacity",
              BAD_CAST "OPACITY" );
  xmlNewProp( curNode,
              BAD_CAST "lineColor",
              BAD_CAST "YELLOW" );
  xmlNewProp( curNode,
              BAD_CAST "id",
              BAD_CAST DefaultNumber );
  curNode = xmlNewChild( curNode, NULL,
                         BAD_CAST "spatialCoordinateCollection", NULL );
  
  // Now that the collection is setup, we iterate through our multimap to
  // populate the list of spatial coordinates
  count = 0;
  lastValue = -1; // worried about this
  for( typename ContourContainerType::const_iterator itr = m_Contours.begin();
       itr != m_Contours.end();
       ++itr )
    {

    // Iterate through each vector in the multimap for the contour points
    unsigned int vectorCount = 0;
    typename ContourPointVectorType::const_iterator vitr;
    for( vitr = itr->second.begin();
         vitr != itr->second.end();
         ++vitr )
      {
      std::string coordinateId = boost::lexical_cast<std::string>(count);
      std::string coordinateIndex = 
        boost::lexical_cast<std::string>( vectorCount );
      ReferenceIndexValueType x = (*vitr)[0];
      ReferenceIndexValueType y = (*vitr)[1];
      ReferenceIndexValueType z = (*vitr)[2];
      std::string strX = boost::lexical_cast<std::string>( x );
      std::string strY = boost::lexical_cast<std::string>( y );
      std::string strZ = boost::lexical_cast<std::string>( z );

      xmlNodePtr coordinateNode = xmlNewChild( curNode,
                                               NULL,
                                               BAD_CAST "SpatialCoordinate",
                                               NULL );
      xmlNewProp( coordinateNode,
                  BAD_CAST "xsi:type",
                  BAD_CAST "TwoDimensionSpacialCoordinate" );
      xmlNewProp( coordinateNode,
                  BAD_CAST "y",
                  BAD_CAST strY.c_str() );
      xmlNewProp( coordinateNode,
                  BAD_CAST "x",
                  BAD_CAST strX.c_str() );
      xmlNewProp( coordinateNode,
                  BAD_CAST "referencedFrameNumber",
                  BAD_CAST strZ.c_str() );
      xmlNewProp( coordinateNode,
                  BAD_CAST "imageReferenceUID",
                  BAD_CAST m_SOPInstanceUIDs.at(z).c_str() );
      xmlNewProp( coordinateNode,
                  BAD_CAST "coordinateIndex",
                  BAD_CAST coordinateIndex.c_str() );
      xmlNewProp( coordinateNode,
                  BAD_CAST "id",
                  BAD_CAST coordinateId.c_str() );

      ++vectorCount;
      ++count;
      }

    }
  
  // Dump the resultant xml to a buffer and into m_Output
  xmlDocDumpFormatMemory(doc, &xmlBuffer, &xmlBufferSize, 1);
  char* docBuffer = (char*) xmlBuffer;
  m_Output = std::string( docBuffer );

  // Free libxml2 memory
  xmlFree(xmlBuffer);
  xmlFreeDoc(doc);
}

//----------------------------------------------------------------------------
template <class TInputImage, class TReferenceImage>
std::string
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::GetCurrentTimeInAIMFormat() const
{
  std::stringstream oss;
  time_t rawtime;
  time( &rawtime );
  struct tm * timeinfo = localtime ( &rawtime );
  oss << 1900+timeinfo->tm_year << "-" 
      << std::setw(2) << setfill('0') << timeinfo->tm_mon << "-"
      << std::setw(2) << setfill('0') << timeinfo->tm_mday << "T"
      << std::setw(2) << setfill('0') << timeinfo->tm_hour << ":"
      << std::setw(2) << setfill('0') << timeinfo->tm_min << ":"
      << std::setw(2) << setfill('0') << timeinfo->tm_sec;
  return oss.str();
}

//----------------------------------------------------------------------------
template <class TInputImage, class TReferenceImage>
void
ImageToAIMXMLFilter<TInputImage,TReferenceImage>
::Update()
{
  this->PopulateContourContainer();
  this->GenerateXMLFromContours();
}

} // end namespace itk

#endif
