//		CTree.h

#ifndef CTREE_H
#define CTREE_H

#include <vector>
#include <velesquant/models/utility.h>

namespace velesquant {

// using namespace std; // REMOVED by refactor
typedef std::vector<double> Vdoub;
typedef std::vector<std::vector<double>> Mdoub;

enum tType { recomb = 0, nonrecomb };
enum exStyle { American = 0, European, Bermudan };

/**
 * @class CTree
 * @brief Binomial and Trinomial Tree Implementation.
 *
 * Used for pricing options using tree methods (Jarrow-Rudd,
 * Cox-Ross-Rubinstein, Trigeorgis).
 */
class CTree {
public:
  CTree();
  CTree(double S, Vdoub T, Vdoub F, Vdoub IV, Vdoub r = Vdoub(),
        Vdoub q = Vdoub());
  ~CTree();
  /**
   * @brief Calculates option price using Binomial Tree.
   *
   * @param strike Option strike.
   * @param Maturity Option maturity.
   * @param Nnodes Number of time steps.
   * @param style Exercise style (American, European, Bermudan).
   * @param pay Call or Put.
   * @param tree Tree type (recombining or non-recombining).
   * @return Option Price.
   */
  double calculateBinomial(double strike, double Maturity, int Nnodes,
                           exStyle style, OptionType pay, tType tree);
  double calculateBinomial(double strike, double Maturity, int Nnodes,
                           exStyle style, OptionType pay, tType tree,
                           Vdoub schedule);

  /**
   * @brief Calculates option price using Trinomial Tree.
   */
  double calculateTrinomial(double strike, double Maturity, int Nnodes,
                            exStyle style, OptionType pay, tType tree);

private:
  double spot_;

  Vdoub T_, F_, IV_, r_, q_;
  // Mdoub	IV_;

  double interp_vol(std::vector<double> &T, std::vector<double> &F, double t1,
                    double t2) const;
  double interpolate(std::vector<double> &T, std::vector<double> &F,
                     double t) const;
};

} // namespace velesquant
#endif