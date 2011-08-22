/*=========================================================================
  Author: $Author: krm15 $  // Author of last commit
  Version: $Rev: 585 $  // Revision of last commit
  Date: $Date: 2009-08-20 21:25:19 -0400 (Thu, 20 Aug 2009) $  // Date of last commit
=========================================================================*/

/*=========================================================================
 Authors: The GoFigure Dev. Team.
 at Megason Lab, Systems biology, Harvard Medical school, 2009

 Copyright (c) 2009, President and Fellows of Harvard College.
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met:

 Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 Redistributions in binary form must reproduce the above copyright notice,
 this list of conditions and the following disclaimer in the documentation
 and/or other materials provided with the distribution.
 Neither the name of the  President and Fellows of Harvard College
 nor the names of its contributors may be used to endorse or promote
 products derived from this software without specific prior written
 permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS
 BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,
 OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT
 OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
 OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/
#if defined(_MSC_VER)
#pragma warning ( disable : 4786 )
#endif

#ifdef __BORLANDC__
#define ITK_LEAN_AND_MEAN
#endif

// Software Guide : BeginLatex
//
// This example illustrates the use of the
// \doxygen{HoughTransform2DSpheresImageFilter} to find circles in a
// N-dimensional image.
//
// First, we include the header files of the filter.
//
// Software Guide : EndLatex


// Software Guide : BeginCodeSnippet
#include "itkHoughTransformRadialVotingImageFilter.h"
// Software Guide : EndCodeSnippet

#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkThresholdImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"
#include <itkGradientMagnitudeImageFilter.h>
#include <itkDiscreteGaussianImageFilter.h>
#include <list>
#include "itkCastImageFilter.h"
#include "vnl/vnl_math.h"

#include "time.h"

int main( int argc, char *argv[] )
{
  if( argc < 15 )
    {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage: " << argv[0] << std::endl;
    std::cerr << " inputImage " << std::endl;
    std::cerr << " accumulatorImage" << std::endl;
    std::cerr << " numberOfSpheres " << std::endl;
    std::cerr << " radius Min " << std::endl;
    std::cerr << " radius Max " << std::endl;
    std::cerr << " SigmaGradient (default = 1) " << std::endl;
    std::cerr << " variance of the accumulator blurring (default = 1) " << std::endl;
    std::cerr << " radius ratio of the disk to remove from the accumulator (default = 1) "<< std::endl;
	  std::cerr << " voting radius ratio (default = 0.5) "<< std::endl;
		std::cerr << " input threshold "<< std::endl;
    std::cerr << " output threshold "<< std::endl;
    std::cerr << " gradient threshold "<< std::endl;
	  std::cerr << " number of threads "<< std::endl;
	  std::cerr << " sampling ratio "<< std::endl;
		return 1;
    }

	clock_t beginT, endT;
  beginT = clock();

  const    unsigned int    Dimension = 3;
  typedef  unsigned char   InputPixelType;
  typedef  float           InternalPixelType;
  typedef  unsigned char   OutputPixelType;

  typedef itk::Image< InputPixelType, Dimension >  InputImageType;
  typedef itk::Image< InternalPixelType, Dimension >    InternalImageType;
  typedef  itk::Image< OutputPixelType, Dimension >  OutputImageType;

  InputImageType::SpacingType spacing;

  typedef  itk::ImageFileReader< InputImageType > ReaderType;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( argv[1] );
  try
    {
    reader->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }
  InputImageType::Pointer localImage = reader->GetOutput();
  spacing = localImage->GetSpacing();
  std::cout << "Computing Hough Map" << std::endl;

  typedef itk::HoughTransformRadialVotingImageFilter< InputImageType,
               InternalImageType > HoughTransformFilterType;
  HoughTransformFilterType::Pointer houghFilter = HoughTransformFilterType::New();
  houghFilter->SetInput( reader->GetOutput() );
  houghFilter->SetNumberOfSpheres( atoi(argv[3]) );
  houghFilter->SetMinimumRadius(   atof(argv[4]) );
  houghFilter->SetMaximumRadius(   atof(argv[5]) );
  houghFilter->SetSigmaGradient( atof(argv[6]) );
  houghFilter->SetVariance( atof(argv[7]) );
  houghFilter->SetSphereRadiusRatio( atof(argv[8]) );
  houghFilter->SetVotingRadiusRatio( atof(argv[9]) );
  houghFilter->SetThreshold( atof(argv[10]) );
  houghFilter->SetOutputThreshold( atof(argv[11]) );
  houghFilter->SetGradientThreshold( atof(argv[12]) );
  houghFilter->SetNbOfThreads( atoi(argv[13]) );
  houghFilter->SetSamplingRatio( atof(argv[14]) );
  houghFilter->Update();

  InternalImageType::Pointer localAccumulator = houghFilter->GetOutput();

  typedef  itk::ImageFileWriter< InternalImageType  > InputWriterType;
  InputWriterType::Pointer writer = InputWriterType::New();

  writer->SetFileName( argv[2] );
  writer->SetInput( localAccumulator );
  writer->SetUseCompression( true );

  try
    {
    writer->Update();
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << "Exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    }

	endT = clock();
	std::cout << ( ( endT - beginT )/CLOCKS_PER_SEC ) << std::endl;

  return 0;
}
