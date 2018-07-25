/******************************************************************************
 * Authors:  Laurent Kneip & Paul Furgale                                     *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

//Note: has been derived from ROS
#include <iostream>

template<typename P>
opengv::sac::Ransac<P>::Ransac(
    int maxIterations, double threshold, double probability) :
    SampleConsensus<P>(maxIterations, threshold, probability)
{}

template<typename P>
opengv::sac::Ransac<P>::~Ransac(){}


template<typename PROBLEM_T>
bool
opengv::sac::Ransac<PROBLEM_T>::computeModel(
    int debug_verbosity_level)
{
  std::cout<<"Inside Compute model"<<std::endl;
  std::cout<<"Type defs"<<std::endl;
  typedef PROBLEM_T problem_t;
  typedef typename problem_t::model_t model_t;

  std::cout<<"var def"<<std::endl;
  iterations_ = 0;
  int n_best_inliers_count = -INT_MAX;
  double k = 1.0;

  std::vector<int> selection;
  model_t model_coefficients;

  int n_inliers_count = 0;
  unsigned skipped_count = 0;
  // supress infinite loops by just allowing 10 x maximum allowed iterations for
  // invalid model parameters!
  const unsigned max_skip = max_iterations_ * 10;

  // Iterate
  std::cout<<"While loop!"<<std::endl;
  while( iterations_ < k && skipped_count < max_skip )
  {
    // Get X samples which satisfy the model criteria
    std::cout<<"While] sample"<<std::endl;
    sac_model_->getSamples( iterations_, selection );

    std::cout<<"While] empty check"<<std::endl;
    if(selection.empty()) 
    {
      fprintf(stderr,
          "[sm::RandomSampleConsensus::computeModel] No samples could be selected!\n");
      break;
    }

    std::cout<<"While] coeff compute"<<std::endl;
    // Search for inliers in the point cloud for the current plane model M
    if(!sac_model_->computeModelCoefficients( selection, model_coefficients ))
    {
      //++iterations_;
      ++ skipped_count;
      continue;
    }

    // Select the inliers that are within threshold_ from the model
    //sac_model_->selectWithinDistance( model_coefficients, threshold_, inliers );
    //if(inliers.empty() && k > 1.0)
    //  continue;

    std::cout<<"While] count dist"<<std::endl;
    n_inliers_count = sac_model_->countWithinDistance(
        model_coefficients, threshold_ );

    // Better match ?
    std::cout<<"While] better match"<<std::endl;
    if(n_inliers_count > n_best_inliers_count)
    {
      std::cout<<"While] inliner"<<std::endl;
      n_best_inliers_count = n_inliers_count;

      // Save the current model/inlier/coefficients selection as being the best so far
      std::cout<<"While] model check"<<std::endl;
      model_              = selection;
      model_coefficients_ = model_coefficients;

      // Compute the k parameter (k=log(z)/log(1-w^n))
      std::cout<<"While] static cast"<<std::endl;
      double w = static_cast<double> (n_best_inliers_count) /
          static_cast<double> (sac_model_->getIndices()->size());
      std::cout<<"While] pow"<<std::endl;
      double p_no_outliers = 1.0 - pow(w, static_cast<double> (selection.size()));
      std::cout<<"While] numeric limits epsilon"<<std::endl;
      p_no_outliers =
          (std::max) (std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by -Inf
      p_no_outliers =
          (std::min) (1.0 - std::numeric_limits<double>::epsilon(), p_no_outliers);
          // Avoid division by 0.
      std::cout<<"While] log"<<std::endl;
      k = log(1.0 - probability_) / log(p_no_outliers);
    }

    ++iterations_;
    std::cout<<"While] verbose log"<<std::endl;
    if(debug_verbosity_level > 1)
      fprintf(stdout,
          "[sm::RandomSampleConsensus::computeModel] Trial %d out of %f: %d inliers (best is: %d so far).\n",
          iterations_, k, n_inliers_count, n_best_inliers_count );
    if(iterations_ > max_iterations_)
    {
      std::cout<<"While] non verbose log"<<std::endl;
      if(debug_verbosity_level > 0)
        fprintf(stdout,
            "[sm::RandomSampleConsensus::computeModel] RANSAC reached the maximum number of trials.\n");
      break;
    }
  }
  std::cout<<"fprintf"<<std::endl;

  if(debug_verbosity_level > 0)
    fprintf(stdout,
        "[sm::RandomSampleConsensus::computeModel] Model: %zu size, %d inliers.\n",
        model_.size(), n_best_inliers_count );

  std::cout<<"model is empty check"<<std::endl;
  if(model_.empty())
  {
    inliers_.clear();
    return (false);
  }

  // Get the set of inliers that correspond to the best model found so far
  std::cout<<"select within distance"<<std::endl;
  sac_model_->selectWithinDistance( model_coefficients_, threshold_, inliers_ );

  return (true);
}
