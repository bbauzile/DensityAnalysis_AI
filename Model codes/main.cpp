// Initially created by Alessio Andronico (alessio.andronico@pasteur.fr) (Andronico et a. 2019)
// adjusted to investigate the impact of reducing the density of palmiped and computed the Re By B. Bauzile and B. Durand

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <iterator>
#include <unordered_set>
#include <unordered_map>
#include <random>
#include <functional>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include <climits>
#include <algorithm>
#include <ctime>
#include <cassert>
#include <iomanip>
using namespace std;


// Constants
const double PI = 3.14159265358979323846;
const double TWOPI = 2.0 * PI;
const int JAN5TH = 38; // First preventive culling (from data)
const int JAN9TH = 42; // Beginning of serious preventive culling
const int FEB1ST = 65; // Beginning of second epidemic wave
const int FEB2ND = 66; // Beginning of systematic preventive culling


// -----------------------------------------------------------------------------
// The parent of all models (abstract base class)
// -----------------------------------------------------------------------------
struct Model {
  size_t nParams; // Total number of parameters in the model
  vector<double> params; // Parameters

  // Virtual methods (must be implemented in subclasses)
  virtual void setParameter(size_t param_i, double value) = 0;
  virtual double computeLogLik() = 0;
  virtual string toString() const = 0;

  // Output
  friend ostream& operator<<(ostream &output, const Model &model) {
    output << model.toString();
    return output;
  }
};


// -----------------------------------------------------------------------------
// MCMC class
// -----------------------------------------------------------------------------
struct MCMC {
  Model *model; // Pointer to model used
  string outSuffix; // Suffix used for output file (file will be named "output" + outSuffix + "dat")
  size_t nIterations;
  size_t thinningFactor;
  double logLik;
  vector<double> lowerLimits;
  vector<double> upperLimits;
  vector<double> priorMeans;
  vector<double> priorStds;
  vector<double> rwSteps; // Random walk steps
  vector<double> acceptanceRates;
  vector<size_t> proposalType; // 0: lognormal proposal, 1: normal proposal
  vector<size_t> priorType; // 0: uniform, 1: gaussian
  mt19937_64 gen; // Mersenne Twister 19937 generator
  normal_distribution<double> normDist;
  uniform_real_distribution<double> unifDist;

  MCMC(Model *model, size_t seed) {
    // Initialize model parameters
    this->model = model;
    outSuffix = "";
    nIterations = 1000;
    thinningFactor = 5;
    logLik = -1.0;
    lowerLimits = vector<double>(model->nParams, -1e20);
    upperLimits = vector<double>(model->nParams, 1e20);
    priorMeans = vector<double>(model->nParams, 0.0);
    priorStds = vector<double>(model->nParams, 1e20);
    rwSteps = vector<double>(model->nParams, 0.5);
    acceptanceRates = vector<double>(model->nParams, 0.0);
    proposalType = vector<size_t>(model->nParams, 0);
    priorType = vector<size_t>(model->nParams, 0);
    // Initialize random number generator seed
    gen.seed(seed);
    // Normal distribution (mean = 0, std = 1)
    normDist = normal_distribution<double>(0.0, 1.0);
    // Uniform distribution in ]0, 1[
    unifDist = uniform_real_distribution<double>(nextafter(0.0, DBL_MAX), 1.0);
  }

  void setBounds(size_t param_i, double lower, double upper) {
    assert(param_i >= 0);
    assert(param_i < model->nParams);
    lowerLimits[param_i] = std::min(lower, upper);
    upperLimits[param_i] = std::max(lower, upper);
    // If current value of parameter is out of bounds, reset it to mid point
    // (otherwise it will never be updated)
    if (model->params[param_i] < lowerLimits[param_i] ||
        model->params[param_i] > upperLimits[param_i]) {
      model->params[param_i] = 0.5 * (lowerLimits[param_i] +
        upperLimits[param_i]);
    }
  }

  double updateNormalProposalUniformPrior(size_t paramNumber) {
    double oldValue = model->params[paramNumber];
    double newValue = oldValue + rwSteps[paramNumber] * normDist(gen);
    double logProposal = 0.0;

    if (newValue < lowerLimits[paramNumber]) return 0.0;
    else if (newValue > upperLimits[paramNumber]) return 0.0;

    model->setParameter(paramNumber, newValue);
    double newLogLik = model->computeLogLik();
    double Q = newLogLik - logLik + logProposal;

    if (log(unifDist(gen)) < Q) {
      logLik = newLogLik;
      return 1.0;
    } else {
      model->setParameter(paramNumber, oldValue);
      return 0.0;
    }
  }

  double updateLogNormalProposalUniformPrior(size_t paramNumber) {
    double oldValue = model->params[paramNumber];
    double exponent = rwSteps[paramNumber] * normDist(gen);
    double newValue = oldValue * exp(exponent);
    double logProposal = exponent; // log(newValue) - log(oldValue);

    if (newValue < lowerLimits[paramNumber]) return 0.0;
    else if (newValue > upperLimits[paramNumber]) return 0.0;

    model->setParameter(paramNumber, newValue);
    double newLogLik = model->computeLogLik();
    double Q = newLogLik - logLik + logProposal;

    if (log(unifDist(gen)) < Q) {
      logLik = newLogLik;
      return 1.0;
    } else {
      model->setParameter(paramNumber, oldValue);
      return 0.0;
    }
  }

  double updateNormalProposalNormalPrior(size_t paramNumber) {
    double oldValue = model->params[paramNumber];
    double newValue = oldValue + rwSteps[paramNumber] * normDist(gen);
    double logProposal = 0.0;
    double logPriors = oldValue - newValue;
    logPriors *= ((oldValue + newValue) * 0.5 - priorMeans[paramNumber]);
    logPriors /= (priorStds[paramNumber] * priorStds[paramNumber]);

    model->setParameter(paramNumber, newValue);
    double newLogLik = model->computeLogLik();
    double Q = newLogLik - logLik + logPriors + logProposal;

    if (log(unifDist(gen)) < Q) {
      logLik = newLogLik;
      return 1.0;
    } else {
      model->setParameter(paramNumber, oldValue);
      return 0.0;
    }
  }

  double updateLogNormalProposalNormalPrior(size_t paramNumber) {
    double oldValue = model->params[paramNumber];
    double exponent = rwSteps[paramNumber] * normDist(gen);
    double newValue = oldValue * exp(exponent);
    double logProposal = exponent; // log(newValue) - log(oldValue);
    double logPriors = oldValue - newValue;
    logPriors *= ((oldValue + newValue) * 0.5 - priorMeans[paramNumber]);
    logPriors /= (priorStds[paramNumber] * priorStds[paramNumber]);

    model->setParameter(paramNumber, newValue);
    double newLogLik = model->computeLogLik();
    double Q = newLogLik - logLik + logPriors + logProposal;

    if (log(unifDist(gen)) < Q) {
      logLik = newLogLik;
      return 1.0;
    } else {
      model->setParameter(paramNumber, oldValue);
      return 0.0;
    }
  }

  inline double updateParameter(size_t iParam) {
    if (proposalType[iParam] == 1) {
      if (priorType[iParam] == 1) {
        return updateNormalProposalNormalPrior(iParam);
      } else {
        return updateNormalProposalUniformPrior(iParam);
      }
    } else {
      if (priorType[iParam] == 1) {
        return updateLogNormalProposalNormalPrior(iParam);
      } else {
        return updateLogNormalProposalUniformPrior(iParam);
      }
    }
  }

  void run() {
    vector<double> nMovesAccepted(model->nParams, 0.0);
    vector<double> nMovesProposed(model->nParams, 0.0);
    ofstream outputFile("output" + outSuffix + ".dat");

    clock_t tBeg = clock();
    logLik = model->computeLogLik();
    cout << setprecision(4) << "LL = " << logLik << endl;
    for (size_t iter = 0; iter < nIterations; ++iter) {
      // cout << iter << "/" << flush;
      for (size_t thin = 0; thin < thinningFactor; ++thin) {
        for (size_t param = 0; param < model->nParams; ++param) {
          if (rwSteps[param] > 0.0) {
            nMovesAccepted[param] += updateParameter(param);
            nMovesProposed[param] += 1.0;
            acceptanceRates[param] =
              nMovesAccepted[param] / nMovesProposed[param];
          }
        }
      }
      // Output
      outputFile << logLik << " ";
      outputFile << *model << " ";
      for (size_t curr = 0; curr < model->nParams - 1; ++curr) {
        outputFile << acceptanceRates[curr] << " ";
      }
      outputFile << acceptanceRates[model->nParams - 1] << endl;
    }
    outputFile.close();
    clock_t tEnd = clock();
    double elapsed_secs = double(tEnd - tBeg) / CLOCKS_PER_SEC;
    // cout << endl;
    cout << "LL = " << logLik << endl;
    cout << "MCMC took: " << elapsed_secs << "s" << endl;
    cout << endl;
    cout << "----- Parameters, rwSteps, and Acceptance Rates ----" << endl;
    for (size_t iParam = 0; iParam < model->nParams; ++iParam) {
      cout << iParam << ": " << fixed << setw(10) << setprecision(7) <<
        model->params[iParam] << ", " << rwSteps[iParam] << ", " <<
        acceptanceRates[iParam] << endl;
    }
    cout << "----------------------------------------------------" << endl;
  }
};


// -----------------------------------------------------------------------------
// Cell class and related functions
// -----------------------------------------------------------------------------
struct Cell {
  string ID;
  string type;
  string department;
  string INSEE;
  double x, y; // Cordinates
  int dateInfection;
  int dateInfectious;
  int dateDetection;
  int dateCulling;
  int dateEntrySZ;
  bool isSymptomatic;
  bool isPreventivelyCulled;
  bool inZRP;
  bool isPAG;
  int lastDayExposure;

  // Constructors
  Cell() { resetInfection(); }

  Cell(string ID, string type, string department, string INSEE,
    double x, double y, int dateInfection, int dateInfectious,
    int dateDetection, int dateCulling, int dateEntrySZ,
    bool isSymptomatic, bool isPreventivelyCulled,
    bool inZRP, bool isPAG):
    ID(ID), type(type), department(department), INSEE(INSEE), x(x), y(y),
    dateInfection(dateInfection), dateInfectious(dateInfectious),
    dateDetection(dateDetection), dateCulling(dateCulling),
    dateEntrySZ(dateEntrySZ), isSymptomatic(isSymptomatic),
    isPreventivelyCulled(isPreventivelyCulled),
    inZRP(inZRP), isPAG(isPAG) { }

  // Copy constructor
  Cell(const Cell& otherCell) { copy(otherCell); }

  // Assignment operator
  Cell& operator=(const Cell& otherCell) {
    if (this != &otherCell) {
      copy(otherCell);
    }
    return *this;
  }

  // Helper function to copy all cell fields
  void copy(const Cell& otherCell) {
    ID = otherCell.ID;
    type = otherCell.type;
    department = otherCell.department;
    INSEE = otherCell.INSEE;
    x = otherCell.x;
    y = otherCell.y;
    dateInfection = otherCell.dateInfection;
    dateInfectious = otherCell.dateInfectious;
    dateDetection = otherCell.dateDetection;
    dateCulling = otherCell.dateCulling;
    dateEntrySZ = otherCell.dateEntrySZ;
    isSymptomatic = otherCell.isSymptomatic;
    isPreventivelyCulled = otherCell.isPreventivelyCulled;
    inZRP = otherCell.inZRP;
    isPAG = otherCell.isPAG;
    lastDayExposure = otherCell.lastDayExposure;
  }

  void resetInfection() {
    dateInfection = INT_MAX;
    dateInfectious = INT_MAX;
    dateDetection = INT_MAX;
    dateCulling = INT_MAX;
    dateEntrySZ = INT_MAX;
    isSymptomatic = false;
    isPreventivelyCulled = false;
    lastDayExposure = INT_MAX;
  }

  // Output
  friend ostream& operator<<(ostream& output, const Cell& cell) {
    int dateDetection = cell.dateDetection;
    int dateInfection = cell.dateInfection;
    int dateInfectious = cell.dateInfectious;
    dateDetection = (dateDetection == INT_MAX) ? -1 : dateDetection;
    dateInfection = (dateInfection == INT_MAX) ? -1 : dateInfection;
    dateInfectious = (dateInfectious == INT_MAX) ? -1 : dateInfectious;
    int dateCulling = cell.dateCulling;
    dateCulling = (dateCulling == INT_MAX) ? -1 : dateCulling;
    output << cell.ID << " " << cell.department << " " << cell.INSEE << " "
      << dateDetection << " " << dateCulling << " "
      << int(cell.isSymptomatic) << " " << int(cell.isPreventivelyCulled)
      << " " << dateInfection << " " << dateInfectious;
    return output;
  }
};


// Euclidean distance between two cells
double distance(const Cell& cellA, const Cell& cellB) {
  double dummyX = cellA.x - cellB.x;
  double dummyY = cellA.y - cellB.y;
  return sqrt(dummyX * dummyX + dummyY * dummyY);
}


// -----------------------------------------------------------------------------
// Space-time survival models similar to:
//   P. Walker et al., J. R. Soc. 9, 1836-1845 (2012)
// -----------------------------------------------------------------------------
struct SurvivalModel: public Model {
  bool preemptive_culling;
  int minDateInfection;
  int minDateDetection;
  int minDateInfectious; // First day an infectious farm exists
  int maxDate; // Max date of culling or detection.
  unordered_map<string, int> latentPeriod; // Infection -> infectiousness
  unordered_map<string, int> incubationPeriod; // Infection -> symptoms
  int detectionDelay; // Delay between symptom onset and detection
  int delayAsymptomatic; // Delay infection -> detection (or culling) for asymptomatic
  int cullingTargets; // 0 = all birds, 1 = palmipedes only, 2 = PAGs only
  double cullingRadius0; // Culling distance around infection (all birds)
  double cullingRadius1; // Culling distance for PAGs
  double surveillanceRadius; // Surveillance radius zone around infection
  double cutoffDistance;
  size_t nCells; // Total number of cells
  vector<Cell> cells; // Container of all cells
  unordered_set<size_t> infectedCellsIDs; // Indices of infected cells
  unordered_set<size_t> survivedCellsIDs; // Indices of cells that escaped infection
  unordered_map<int, unordered_set<size_t>> infectiousCellsAt; // Indices of cells infectious at a given time
  vector<vector<double>> distanceMatrix; // Distance matrix between cells

  size_t nSwitchPoints;
  vector<int> switchPoints;

  // To keep track of of cells within certain distances
  unordered_map<size_t, unordered_set<size_t>> noPAGCellsWithinCullingRadius0;
  unordered_map<size_t, unordered_set<size_t>> cellsWithinCullingRadius0;
  unordered_map<size_t, unordered_set<size_t>> PAGCellsWithinCullingRadius1;
  unordered_map<size_t, unordered_set<size_t>> cellsWithinSurveillanceRadius;
  unordered_map<size_t, unordered_set<size_t>> cellsWithinCutoffDistance;
  vector<size_t> nWithinCutoffDistance;

  unordered_set<string> departments;

  // Random numbers and distributions
  double probaPassiveIP;
  double probaPrevCulledIP;
  double probaActiveIP;
  mt19937_64 gen;
  uniform_real_distribution<double> unifDist;
  gamma_distribution<double> gammaCull;
  gamma_distribution<double> gammaDet;
  gamma_distribution<double> gammaPrevCull;
  discrete_distribution<size_t> choiceIP;
  discrete_distribution<size_t> choiceIP_NoPreventive;

  SurvivalModel(string dataFile,
      int latentPeriod_Palmi, int latentPeriod_Galli,
      int incubationPeriod_Palmi, int incubationPeriod_Galli,
      int detectionDelay, int delayAsymptomatic,
      double cutoffDistance, const vector<int>& switchPoints) {
    assert(incubationPeriod_Palmi >= latentPeriod_Palmi);
    assert(incubationPeriod_Galli >= latentPeriod_Galli);
    this->latentPeriod["palmipedes"] = latentPeriod_Palmi;
    this->latentPeriod["galliformes"] = latentPeriod_Galli;
    this->incubationPeriod["palmipedes"] = incubationPeriod_Palmi;
    this->incubationPeriod["galliformes"] = incubationPeriod_Galli;
    this->detectionDelay = detectionDelay;
    this->delayAsymptomatic = delayAsymptomatic;
    this->switchPoints = switchPoints;
    nSwitchPoints = switchPoints.size();

    // From data
    probaPassiveIP = 0.44;
    probaPrevCulledIP = 0.30;
    probaActiveIP = 0.26;
    choiceIP = discrete_distribution<size_t>{
      probaPassiveIP, probaActiveIP, probaPrevCulledIP
    };
    choiceIP_NoPreventive = discrete_distribution<size_t>{
      probaPassiveIP, probaActiveIP
    };
    // Delay detection -> culling
    gammaCull = gamma_distribution<double>(2.4376128, 1.0 / 0.4698795);
    // Delay culling -> detection (for preventively culled (asymptomatic) infections)
    gammaDet = gamma_distribution<double>(2.50888690, 1.0 / 0.64574747);
    // Delay preventive culling (from detection date)
    gammaPrevCull = gamma_distribution<double>(1.67931803, 1.0 / 0.21506603);

    // Initialize distributions
    gen.seed(0);
    unifDist = uniform_real_distribution<double>(0.0, nextafter(1.0, DBL_MAX));

    loadData(dataFile);

    cullingTargets = 0;
    initializeDistanceMatrix();
    setCullingRadius0(1.0);
    setCullingRadius1(3.0);
    setSurveillanceRadius(10.0);
    setCutoffDistance(cutoffDistance);

    assignHistory(cells, infectiousCellsAt, infectedCellsIDs, survivedCellsIDs);
    for (int switchPoint : switchPoints) {
      assert(switchPoint > minDateInfection && switchPoint < maxDate);
    }
    verifyData();

    nParams = 5 + nSwitchPoints + 1;
    params = vector<double>(nParams, 1.0);

    preemptive_culling = true;
  }

  void loadData(string dataFile) {
    // Assumptions about file structure and data contained in it:
    // 1) 14 columns (space separated): ID, department, INSEE, type, x, y,
    //      dateDetection, dateCulling, Symptomatic, PreventiveCulling,
    //      delayEmpty, coop, isPAG, inZRP
    // 2) dateDetection is >= 0 for infected cells
    // 3) Missing dates are = -1
    // 4) type is either "palmipedes" or "galliformes"
    ifstream iFile(dataFile);

    // Check that file is there
    if (!iFile.is_open()) {
      cout << "ERROR: file " << dataFile << " does not exist!" << endl;
      exit(1);
    }

    minDateInfection = INT_MAX;
    minDateDetection = INT_MAX;

    // Read data
    string ID;
    string department;
    string INSEE;
    string type;
    double x, y;
    int dateInfection;
    int dateInfectious;
    int dateDetection;
    int dateCulling;
    int symptomatic;
    int prevCulling;
    int delayEmpty;
    string coop;
    int dateEntrySZ = INT_MAX;
    bool isSymptomatic;
    bool isPreventivelyCulled;
    int isPAG;
    int inZRP;
    size_t curr = 0;
    maxDate = INT_MIN;
    unordered_set<string> typeCheck;
    string line;
    while (!iFile.eof()) {
      getline(iFile, line);
      if (line == "") continue; // Nothing to parse
      istringstream parser(line);
      parser >> ID >> department >> INSEE >> type >> x >> y >> dateDetection
        >> dateCulling >> symptomatic >> prevCulling >> delayEmpty >> coop
        >> isPAG >> inZRP;
      dateDetection = dateDetection == -1 ? INT_MAX : dateDetection;
      dateCulling = dateCulling == -1 ? INT_MAX : dateCulling;
      isSymptomatic = symptomatic == 1;
      isPreventivelyCulled = prevCulling == 1;

      if (dateDetection < INT_MAX) {
        if (isSymptomatic) {
          dateInfection = dateDetection -
            (incubationPeriod[type] + detectionDelay);
        } else {
          if (isPreventivelyCulled) {
            // There are cases in the data where dateDetection < dateCulling
            dateInfection = std::min(dateDetection, dateCulling)
              - delayAsymptomatic;
          } else {
            dateInfection = dateDetection - delayAsymptomatic;
          }
        }
        dateInfectious = dateInfection + latentPeriod[type];
        if (dateInfectious >= dateCulling) {
          cout << "WARNING: foyer " << ID << " will not be infectious" << endl;
          dateInfectious = INT_MAX;
        }
        infectedCellsIDs.insert(curr);
      } else {
        dateInfection = INT_MAX;
        dateInfectious = INT_MAX;
        survivedCellsIDs.insert(curr);
      }

      Cell currCell(ID, type, department, INSEE, x, y, dateInfection,
        dateInfectious, dateDetection, dateCulling, dateEntrySZ,
        isSymptomatic, isPreventivelyCulled, bool(inZRP), bool(isPAG));
      cells.push_back(currCell);

      typeCheck.insert(type);
      departments.insert(department);

      // Keep track of infection and detection dates' range
      if (dateInfection < minDateInfection) minDateInfection = dateInfection;
      if (dateDetection < minDateDetection) minDateDetection = dateDetection;
      if (dateDetection < INT_MAX) {
        maxDate = std::max(maxDate, dateDetection);
      }
      if (dateCulling < INT_MAX) {
        maxDate = std::max(maxDate, dateCulling);
      }
      ++curr;
    }
    iFile.close();
    nCells = cells.size();

    // Check input data
    assert(maxDate >= 0);
    bool check1, check2;
    for (string currType : typeCheck) {
      check1 = currType != "palmipedes";
      check2 = currType != "galliformes";
      if (check1 && check2) {
        cout << "ERROR: got some weird types:" << endl;
        for (string type : typeCheck) {
          cout << type << " ";
        }
        cout << endl;
        cout << "ABORTING NOW!" << endl;
        exit(1);
      }
    }
  }

  void initializeDistanceMatrix() {
    distanceMatrix.resize(nCells, vector<double>(nCells, 0.0));
    for (size_t currI = 0; currI < (nCells - 1); ++currI) {
      for (size_t currJ = currI + 1; currJ < nCells; ++currJ) {
        double dist = distance(cells[currI], cells[currJ]);
        distanceMatrix[currI][currJ] = dist;
        distanceMatrix[currJ][currI] = dist;
      }
    }
  }

  void setCullingRadius0(double radius) {
    assert(radius > 0.0);
    cullingRadius0 = radius;
    cellsWithinCullingRadius0.clear();
    noPAGCellsWithinCullingRadius0.clear();
    for (size_t cellI = 0; cellI < nCells; ++cellI) {
      for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
        if (cellJ != cellI && distanceMatrix[cellI][cellJ] < radius) {
          cellsWithinCullingRadius0[cellI].insert(cellJ);
          if (!cells[cellJ].isPAG) {
            noPAGCellsWithinCullingRadius0[cellI].insert(cellJ);
          }
        }
      }
    }
  }

  void setCullingRadius1(double radius) {
    assert(radius >= cullingRadius0);
    cullingRadius1 = radius;
    PAGCellsWithinCullingRadius1.clear();
    for (size_t cellI = 0; cellI < nCells; ++cellI) {
      for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
        if (cellJ != cellI && cells[cellJ].isPAG
            && distanceMatrix[cellI][cellJ] < radius) {
          PAGCellsWithinCullingRadius1[cellI].insert(cellJ);
        }
      }
    }
  }

  void setSurveillanceRadius(double radius) {
    surveillanceRadius = radius;
    cellsWithinSurveillanceRadius.clear();
    for (size_t cellI = 0; cellI < nCells; ++cellI) {
      for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
        if (cellJ != cellI && distanceMatrix[cellI][cellJ] < radius) {
          cellsWithinSurveillanceRadius[cellI].insert(cellJ);
        }
      }
    }
  }

  void setCutoffDistance(double cutoffDistance) {
    assert(cutoffDistance > 0.0);
    this->cutoffDistance = cutoffDistance;
    cellsWithinCutoffDistance.clear();
    nWithinCutoffDistance.resize(nCells);
    for (size_t cellI = 0; cellI < nCells; ++cellI) {
      for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
        if (cellJ != cellI && distanceMatrix[cellI][cellJ] < cutoffDistance) {
          cellsWithinCutoffDistance[cellI].insert(cellJ);
        }
      }
      nWithinCutoffDistance[cellI] = cellsWithinCutoffDistance[cellI].size();
    }
  }

  void assignHistory(vector<Cell>& vCells,
      unordered_map<int, unordered_set<size_t>>& vInfectiousAt,
      const unordered_set<size_t>& vInfected,
      const unordered_set<size_t>& vSurvived) {
    minDateInfectious = INT_MAX;

    // Compute infectious cells for each date
    vInfectiousAt.clear(); // Make sure map was not already populated
    for (size_t cellI : vInfected) {
      int infectiousI = vCells[cellI].dateInfectious;
      int cullingI = vCells[cellI].dateCulling;
      int stopInfectiousI = std::min(cullingI, maxDate + 1);
      for (int currT = infectiousI; currT < stopInfectiousI; ++currT) {
        vInfectiousAt[currT].insert(cellI);
      }
      if (minDateInfectious > infectiousI) {
        minDateInfectious = infectiousI;
      }
    }

    // Assign dates entry in surveillance zone
    // For each cell keep date of first infection detected within surveillanceRadius
    for (size_t cellI : vInfected) {
      int detectionI = vCells[cellI].dateDetection;
      int currEntryI = vCells[cellI].dateEntrySZ;
      vCells[cellI].dateEntrySZ = std::min(currEntryI, detectionI);
      for (size_t cellJ : cellsWithinSurveillanceRadius[cellI]) {
        int currEntryJ = vCells[cellJ].dateEntrySZ;
        vCells[cellJ].dateEntrySZ = std::min(currEntryJ, detectionI);
      }
    }

    // Assign last day of exposure for all cells
    for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
      int stop = std::min(
        vCells[cellJ].dateInfection, vCells[cellJ].dateCulling
      );
      vCells[cellJ].lastDayExposure = std::min(stop - 1, maxDate);
    }
  }

  void verifyData() {
    // 'assignHistory' must be called before this function
    bool invalidData = false;
    size_t minSize = 1;
    if (latentPeriod["palmipedes"] == 0 || latentPeriod["galliformes"] == 0) {
      minSize = 2;
    }
    for (size_t cellI : infectedCellsIDs) {
      if (cells[cellI].dateInfection == minDateInfection) continue; // Index cells
      int dateInfection = cells[cellI].dateInfection;
      if (infectiousCellsAt[dateInfection].size() < minSize) {
        cout << "WARNING: cell " << cells[cellI].ID
          << " has no possible infectors!" << endl;
        // invalidData = true;
      }
      if (cells[cellI].dateCulling == INT_MAX) {
        cout << "ERROR: cell " << cells[cellI].ID << " has no culling date!"
          << endl;
        invalidData = true;
      }
    }
    if (invalidData) exit(1);
  }

  size_t periodShift(int currDate) {
    size_t index = 0;
    if (nSwitchPoints > 0) {
      auto lower = lower_bound(
        switchPoints.cbegin(), switchPoints.cend(), currDate
      );
      index = distance(switchPoints.cbegin(), lower);
    }
    return index;
  }

  double getPhi(string farmType) {
    // Relative susceptibility of galliformes
    return (farmType == "palmipedes") ? 1.0 : params[0];
  }

  double getPsi(string farmType) {
    // Relative infectivity of galliformes
    return (farmType == "palmipedes") ? 1.0 : params[1];
  }

  double getAlphaSZ(const Cell& cellJ, const Cell& cellI, int day) {
    // Effect of surveillance zone
    if (day >= cellJ.dateEntrySZ || day >= cellI.dateEntrySZ) {
      // Either j or i are in SZ
      return params[2];
    }
    return 1.0;
  }

  double getBetaExt() {
    // External force of infection
    return params[3];
  }

  double getBeta(int currDate, const string& department) {
    // Transmission rate
    if (department != "40") return params[4];
    return params[5 + periodShift(currDate)];
  }

  double foi(size_t cellJ, size_t cellI, int day, const vector<Cell>& vCells) {
    // Force of infection
    if (distanceMatrix[cellI][cellJ] >= cutoffDistance) return 0.0;
    if (day > vCells[cellJ].dateInfection) return 0.0;
    if (day >= vCells[cellJ].dateCulling) return 0.0;
    if (nWithinCutoffDistance[cellI] == 0) return 0.0;
    // Infection possible?
    int infectiousI = vCells[cellI].dateInfectious;
    int cullingI = vCells[cellI].dateCulling;
    if (day >= infectiousI && day < cullingI) {
      double phi = getPhi(vCells[cellJ].type);
      double psi = getPsi(vCells[cellI].type);
      double alphaSZ = getAlphaSZ(vCells[cellJ], vCells[cellI], day);
      double beta = getBeta(day, vCells[cellI].department);
      return beta * alphaSZ * phi * psi / nWithinCutoffDistance[cellI];
    }
    return 0.0;
  }

  double foiExt(size_t cellJ, int day, const vector<Cell>& vCells) {
    // External force of infection
    if (day > vCells[cellJ].dateInfection) return 0.0;
    if (day >= vCells[cellJ].dateCulling) return 0.0;
    double beta = getBetaExt();
    double phi = getPhi(vCells[cellJ].type);
    return beta * phi;
  }

  double getProbaInfection(size_t cellJ, int currDate, vector<Cell>& vCells,
      unordered_map<int, unordered_set<size_t>>& vInfectiousAt) {
    double dummy = 0.0;
    for (size_t cellI : vInfectiousAt[currDate]) {
      if (cellI == cellJ) {
        cout << "ERROR: should never have gotten here (getProbaInfection)!!!"
             << endl;
        exit(1);
      }
      dummy += foi(cellJ, cellI, currDate, vCells);
    }
    dummy += foiExt(cellJ, currDate, vCells);
    double res;
    if (dummy < 1e-20) { // To prevent problems with log
      res = 1e-20;
    } else {
      // If exponent is small, use power series, otherwise full expression
      res = (dummy < 1e-6) ? dummy : 1.0 - exp(-dummy);
    }
    return res;
  }

  // Returns the number of the randomly selected infector from
  // the infectious, or vCells.size() + 1 if the infection is external
  // Must be used before updating vInfectiousAt
  size_t chooseInfector(size_t cellJ, int currDate, vector<Cell>& vCells,
      unordered_map<int, unordered_set<size_t>>& vInfectiousAt) {
    double dummy = 0.0;
    vector<double> all_fois;
    for (size_t cellI : vInfectiousAt[currDate]) {
      if (cellI == cellJ) {
        cout << "ERROR: should never have gotten here (getProbaInfection)!!!"
             << endl;
        exit(1);
      }
      all_fois.push_back(foi(cellJ, cellI, currDate, vCells));
    }
    all_fois.push_back(foiExt(cellJ, currDate, vCells));
    for (auto & foival: all_fois) dummy += foival;
    for (auto & foival: all_fois) foival /= dummy;
    discrete_distribution<size_t> choice(all_fois.begin(), all_fois.end());
    size_t infectorI = choice(gen);
    return (infectorI == all_fois.size() - 1 ? vCells.size() + 1 :
                *(std::next(vInfectiousAt[currDate].begin(), infectorI)));
  }

  double computeLogLik() {
    double loglik = 0.0;

    // Contribution from infections
    for (size_t cellI : infectedCellsIDs) {
      loglik += log(getProbaInfection(cellI, cells[cellI].dateInfection,
        cells, infectiousCellsAt));
    }

    // Contribution from avoided infections
    for (size_t cellI : infectedCellsIDs) {
      double psi = getPsi(cells[cellI].type);
      int cullingI = cells[cellI].dateCulling;
      double den = nWithinCutoffDistance[cellI];
      for (size_t cellJ : cellsWithinCutoffDistance[cellI]) {
        double phi = getPhi(cells[cellJ].type);
        int stop = std::min(cells[cellJ].lastDayExposure, cullingI - 1);
        for (int day = cells[cellI].dateInfectious; day <= stop; ++day) {
          double alphaSZ = getAlphaSZ(cells[cellJ], cells[cellI], day);
          double beta = getBeta(day, cells[cellI].department);
          loglik -= beta * alphaSZ * phi * psi / den;
        }
      }
    }

    // External force of infection
    double beta = getBetaExt();
    for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
      loglik -= beta * (cells[cellJ].lastDayExposure + 1 - minDateInfection)
        * getPhi(cells[cellJ].type);
    }

    return loglik;
  }

  double computeLogLikV2() {
    // Different version of computeLogLik (USED FOR TESTING ONLY)
    double loglik = 0.0;

    // Contribution from infections
    for (size_t cellI : infectedCellsIDs) {
      loglik += log(getProbaInfection(cellI, cells[cellI].dateInfection,
        cells, infectiousCellsAt));
    }

    // Contribution from avoided infections
    for (size_t cellI : infectedCellsIDs) {
      for (size_t cellJ : cellsWithinCutoffDistance[cellI]) {
        for (int day = cells[cellI].dateInfectious;
            day <= cells[cellJ].lastDayExposure; ++day) {
          loglik -= foi(cellJ, cellI, day, cells);
        }
      }
    }

    // External force of infection
    for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
      for (int day = minDateInfection; day <= cells[cellJ].lastDayExposure;
          ++day) {
        loglik -= foiExt(cellJ, day, cells);
      }
    }

    return loglik;
  }

  void rescaleGamma(string gamma_dist, double new_mean) {
    double beta;
    if (gamma_dist == "culling") {
      beta = gammaCull.beta();
      gammaCull = gamma_distribution<double>(new_mean / beta, beta);
    } else if (gamma_dist == "detection") {
      beta = gammaDet.beta();
      gammaDet = gamma_distribution<double>(new_mean / beta, beta);
    } else if (gamma_dist == "preventive_culling") {
      beta = gammaPrevCull.beta();
      gammaPrevCull = gamma_distribution<double>(new_mean / beta, beta);
    } else {
      cout << "WARNING: Skipping rescaling of invalid gamma distribution ("
        << gamma_dist << ")!" << endl;
    }
  }

  void writeFullEpidemic(const vector<Cell>& vCells, ofstream& outputFile) {
    for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
      outputFile << vCells[cellJ] << endl;
    }
  }

  void writeEpidemic(const vector<Cell>& vCells, ofstream& outputFile,
      size_t simId) {
    bool wasInfected;
    bool wasCulled;
    for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
      wasInfected = vCells[cellJ].dateDetection < INT_MAX;
      wasCulled = vCells[cellJ].dateCulling < INT_MAX;
      if (wasInfected || wasCulled) {
        outputFile << simId << " " << vCells[cellJ] << endl;
      }
    }
  }

  void writeTree(const vector<Cell>& vCells,
                 unordered_map<int, unordered_map<size_t, size_t>> & tTree,
                 string fTree, size_t simId)
  {
    ofstream outputFile(fTree, ios::app);
    for (auto & i: tTree) {
        int tps = i.first;
        for (auto & j: i.second)
            outputFile << simId << "," << tps << ","
                       << vCells[j.first].ID << ","
                       << vCells[j.first].type << ","
                       << (j.second < vCells.size() ? vCells[j.second].ID :
                                "external") << ","
                       << (j.second < vCells.size() ? vCells[j.second].type :
                                "external")
                       << endl;
    }
    outputFile.close();
  }

  void copyEpidemicState(int upToDate, vector<Cell>& simCells,
      unordered_set<size_t>& simInfected, unordered_set<size_t>& simSurvived) {
    // NOTE: simCells must be allocated before calling this function
    bool detectedBefore;
    bool culledBefore;
    simInfected.clear();
    simSurvived.clear();
    for (size_t cellJ = 0; cellJ < nCells; ++cellJ) {
      simCells[cellJ] = cells[cellJ];
      detectedBefore = cells[cellJ].dateDetection <= upToDate;
      culledBefore = cells[cellJ].dateCulling <= upToDate;
      if (!(detectedBefore || culledBefore)) {
        simCells[cellJ].resetInfection();
      }
      if (simCells[cellJ].dateInfection < INT_MAX) {
        simInfected.insert(cellJ);
      } else {
        simSurvived.insert(cellJ);
      }
    }
  }

  double probaPreventiveCullingNoPAG0(int date, mt19937_64& gen) {
    // Empirical distribution from data
    if (date < JAN5TH) return 0.0;

    vector<double> probas;
    if (date >= JAN5TH && date < FEB2ND) {
      probas = {0.87012987, 0.01298701, 0.00000000, 0.01298701,
        0.03896104, 0.03896104, 0.00000000, 0.00000000, 0.00000000,
        0.02597403};
    } else {
      probas = {
        0.752475248, 0.009900990, 0.009900990, 0.034653465, 0.064356436,
        0.064356436, 0.004950495, 0.000000000, 0.000000000, 0.059405941
      };
    }
    discrete_distribution<size_t> choice(probas.begin(), probas.end());
    return 0.05 + 0.1 * choice(gen);
  }

  double probaPreventiveCullingPAG1(int date, mt19937_64& gen) {
    // Empirical distribution from data
    if (date < JAN5TH) return 0.0;

    vector<double> probas;
    if (date >= JAN5TH && date < FEB2ND) {
      probas = {
        0.22916667, 0.06250000, 0.14583333, 0.11458333, 0.12500000, 0.11458333,
        0.07291667, 0.08333333, 0.02083333, 0.03125000
      };
    } else {
      probas = {
        0.10294118, 0.04779412, 0.08823529, 0.15441176, 0.21323529, 0.16544118,
        0.08088235, 0.07720588, 0.02573529, 0.04411765
      };
    }
    discrete_distribution<size_t> choice(probas.begin(), probas.end());
    return 0.05 + 0.1 * choice(gen);
  }

  void simulationStep(int currT, vector<Cell>& vCells,
      unordered_set<size_t>& vInfected, unordered_set<size_t>& vSurvived,
      unordered_map<int, unordered_set<size_t>>& vInfectiousAt,
      // BD-20220114
      string fTree,
      unordered_map<int, unordered_map<size_t, size_t>>& tTree,
      unordered_set<size_t> emptied) {

    int detectionI, detectionJ, cullingJ;
    bool notAlreadyCulled;
    double probaInfection;

    unordered_set<size_t> infectedToday;
    for (size_t cellJ : vSurvived) {
      // Cell already culled?
      if (currT > vCells[cellJ].dateCulling) continue;

      // BD-20220114 Cell emptied ?
      if (emptied.find(cellJ) != emptied.end()) continue;

      probaInfection = getProbaInfection(cellJ, currT, vCells, vInfectiousAt);
      if (unifDist(gen) < probaInfection) { // Cell got infected
        // BD-20220114
        if (fTree != "")
            tTree[currT][cellJ] =
                chooseInfector(cellJ, currT, vCells, vInfectiousAt);
        // Assign dates
        vCells[cellJ].dateInfection = currT;
        vCells[cellJ].dateInfectious = currT + latentPeriod[vCells[cellJ].type];
        vCells[cellJ].lastDayExposure = currT - 1;
        bool check_type =
          (vCells[cellJ].type == "galliformes" && cullingTargets == 0) ||
          (vCells[cellJ].type == "palmipedes" &&
            (cullingTargets == 0 || cullingTargets == 1)) ||
          (vCells[cellJ].isPAG && cullingTargets == 2);
        size_t ip_type;
        if (check_type && preemptive_culling) {
          ip_type = choiceIP(gen);
        } else {
          ip_type = choiceIP_NoPreventive(gen);
        }
        if (ip_type == 0) {
          // Passive surveillance
          detectionJ = currT + incubationPeriod[vCells[cellJ].type]
            + detectionDelay;
          vCells[cellJ].isSymptomatic = true;
        } else if (ip_type == 1) {
          // Active surveillance
          detectionJ = currT + delayAsymptomatic;
        } else if (ip_type == 2) {
          // Culled then detected
          cullingJ = currT + delayAsymptomatic;
          vCells[cellJ].dateCulling = cullingJ;
          vCells[cellJ].isPreventivelyCulled = true;
          detectionJ = cullingJ + gammaDet(gen);
        } else {
          cout << "ERROR: should never see this!!" << endl;
          exit(911);
        }
        vCells[cellJ].dateDetection = detectionJ;
        // Reassign culling date (only if not assigned already)
        if (vCells[cellJ].dateCulling == INT_MAX) {
          cullingJ = detectionJ + gammaCull(gen);
          vCells[cellJ].dateCulling = cullingJ;
          vCells[cellJ].lastDayExposure = cullingJ - 1;
        }
        // Assign entry in surveillance zone for all cells around infection
        // (if not already assigned)
        if (vCells[cellJ].dateEntrySZ == INT_MAX) {
          vCells[cellJ].dateEntrySZ = detectionJ;
        }
        for (size_t cellK : cellsWithinSurveillanceRadius[cellJ]) {
          if (vCells[cellK].dateEntrySZ == INT_MAX) {
            vCells[cellK].dateEntrySZ = detectionJ;
          }
        }
        // Add to set of cells infected on this day
        infectedToday.insert(cellJ);
      }
    }
    // Update infection flags
    for (size_t cellI : infectedToday) {
      // Add cell to set of infected cells
      vInfected.insert(cellI);
      // Flag cell as infectious until dateCulling
      for (int tStep = vCells[cellI].dateInfectious;
          tStep < vCells[cellI].dateCulling; ++tStep) {
        vInfectiousAt[tStep].insert(cellI);
      }
      // Remove infected cells from set of survived cells
      vSurvived.erase(cellI);
    }

    // Preventive culling
    if (preemptive_culling) {
      for (size_t cellI : infectedToday) {
        detectionI = vCells[cellI].dateDetection;

        // Culling (no PAG)
        if (cullingTargets == 0 || cullingTargets == 1) {
          for (size_t cellJ : noPAGCellsWithinCullingRadius0[cellI]) {
            if (cullingTargets == 1 && vCells[cellJ].type == "galliformes") {
              continue;
            }
            notAlreadyCulled = vCells[cellJ].dateCulling == INT_MAX;
            if (unifDist(gen) < probaPreventiveCullingNoPAG0(detectionI, gen) &&
                notAlreadyCulled) {
              // Preventive culling
              cullingJ = detectionI + gammaPrevCull(gen);
              vCells[cellJ].dateCulling = cullingJ;
              vCells[cellJ].lastDayExposure = cullingJ - 1;
              vCells[cellJ].isPreventivelyCulled = true;
            }
          }
        }

        // Culling for PAG
        for (size_t cellJ : PAGCellsWithinCullingRadius1[cellI]) {
          notAlreadyCulled = vCells[cellJ].dateCulling == INT_MAX;
          if (unifDist(gen) < probaPreventiveCullingPAG1(detectionI, gen) &&
              notAlreadyCulled) {
            // Preventive culling
            cullingJ = detectionI + gammaPrevCull(gen);
            vCells[cellJ].dateCulling = cullingJ;
            vCells[cellJ].lastDayExposure = cullingJ - 1;
            vCells[cellJ].isPreventivelyCulled = true;
          }
        }
      }
    }
  }

  void simulate(ofstream& outputFile, size_t simID, int fromDate = INT_MIN,
      int toDate = 365,
      //             fTree: name of file to write transmission tree
      //             ePalmi: number of palmi farms to remove to obtain the density threshold of the scenario in the set municipality
      //             eGalli: number of galli farms to remove to obtain the density threshold of the scenario  in the set municipality = 0
      string fTree = "",
      unordered_map<string, size_t> ePalmi = unordered_map<string, size_t>(),
      unordered_map<string, size_t> eGalli = unordered_map<string, size_t>()) {

    // Initialize cells (all data coming <= upToDate are kept)
    int upToDate = (fromDate == INT_MIN) ? minDateInfection : fromDate - 1;
    vector<Cell> vCells(nCells, Cell());
    unordered_set<size_t> vInfected;
    unordered_set<size_t> vSurvived;
    unordered_map<int, unordered_set<size_t>> vInfectiousAt;
    copyEpidemicState(upToDate, vCells, vInfected, vSurvived);
    assignHistory(vCells, vInfectiousAt, vInfected, vSurvived);

    int currT;
    int maxLatent = std::max(latentPeriod["palmipedes"],
      latentPeriod["galliformes"]);
    int maxIncubation = std::max(incubationPeriod["palmipedes"],
      incubationPeriod["galliformes"]);


    unordered_set<size_t> emptied;
    for (auto & i: ePalmi) {
        string insee = i.first;
        size_t nemptied = i.second;
        vector<size_t> farms;
        for (unsigned i = 0 ; i < nCells ; i++)
            if ((vCells[i].INSEE == insee) && (vCells[i].type == "palmipedes"))
                farms.push_back((size_t)i);
        if (farms.size() < nemptied)
            cout << "Number of palmipeds farms is insufficient"
                 << "in the municipality '" << insee << "' (" << farms.size()
                 << ") to be able to remove from it " << nemptied << "\n";
        shuffle(farms.begin(), farms.end(), gen);
        for (unsigned i = 0 ; i < std::min(farms.size(), nemptied) ; i++)
            emptied.insert(farms[i]);
    }
    for (auto & i: eGalli) {
        string insee = i.first;
        size_t nemptied = i.second;
        vector<size_t> farms;
        for (unsigned i = 0 ; i < nCells ; i++)
            if ((vCells[i].INSEE == insee) && (vCells[i].type == "galliformes"))
                farms.push_back((size_t)i);
        if (farms.size() < nemptied)
            cout << "Number of galliforms is insufficient"
                 << "in the municipality '" << insee << "' (" << farms.size()
                 << ") to be able to remove from it " << nemptied << "\n";
        shuffle(farms.begin(), farms.end(), gen);
        for (unsigned i = 0 ; i < std::min(farms.size(), nemptied) ; i++)
            emptied.insert(farms[i]);
    }

    if (fromDate == INT_MIN) {
      // Start from first day with an infectious farm
      currT = minDateInfection
        + std::min(latentPeriod["palmipedes"], latentPeriod["galliformes"]);
    } else {
      // Start from earliest date of infection for a farm detected on day
      // 'fromDate'
      int maxInfectionDelay = std::max(maxIncubation + detectionDelay,
        delayAsymptomatic);
      currT = fromDate - maxInfectionDelay;
    }

    // initialize tTree with the infected farms at t0
    unordered_map<int, unordered_map<size_t, size_t>> tTree;
    for (auto & cellJ: vInfected)
        tTree[currT][cellJ] = vCells.size() + 1;

    while (true) { // Repeat loop until no infectious cells or stopping time
      size_t nInfected = vInfected.size();
      // Time to stop: currT is > stopping time or all cells got infected
      if (currT > toDate || nInfected == nCells) break;

      // Advance simulation
      simulationStep(currT, vCells, vInfected, vSurvived, vInfectiousAt,
                     fTree, tTree, emptied);

      // Check that there are still infectious farms
      size_t daysWithInfectiousFarms = 0;
      for (int thisT = currT + 1; thisT <= (currT + maxLatent); ++thisT) {
        daysWithInfectiousFarms += vInfectiousAt.count(thisT);
      }
      if (daysWithInfectiousFarms == 0) break;

      ++currT;
    }
    writeEpidemic(vCells, outputFile, simID);
    if (fTree != "") writeTree(vCells, tTree, fTree, simID);
  }

  void simulateEpidemics(string mcmc_file, string out_file,
      size_t burn_in, size_t n_params_sets, size_t n_draws,
      int fromDate = (INT_MIN + 1), int toDate = 365, size_t seed = 1987611,
      // BD-20220114 fTree: name of file where to write trees
      // eFarm: file with 3 columns separated by blanks and
      // without column headers: municipality ID, palmi= number
      // of palmi farms to emptied at the start of the simulation,
      // galli: same for galliformes
      string fTree = "", string eFarm = "") {

    if (fromDate > maxDate) {
      cout << "ERROR: fromDate > maxDate!" << endl;
      exit(1);
    }

    // reading eFarm
    unordered_map<string, size_t> ePalmi, eGalli;
    if (eFarm != "") {
        ifstream iFile(eFarm);
        string line = "", insee;
        size_t npalmi = 0, ngalli = 0;
        if (!iFile.is_open()) {
          cout << "ERROR: file " << eFarm << " does not exist!" << endl;
          exit(1);
        }
        while (!iFile.eof()) {
          getline(iFile, line);
          if (line == "") continue; // Nothing to parse
          istringstream parser(line);
          parser >> insee >> npalmi >> ngalli;
          if (npalmi > 0) ePalmi[insee] = npalmi;
          if (ngalli > 0) eGalli[insee] = ngalli;
        }
        iFile.close();

    }

    // Check that MCMC output file is there
    ifstream m_file(mcmc_file);
    if (!m_file.is_open()) {
      cout << "ERROR: file " << mcmc_file << " does not exist!" << endl;
      exit(1);
    }

    // Read in file (leaving out burn_in iterations)
    // NOTE: need to use stringstream pointers to run on cluster (there is
    //  a bug on the gcc version currently installed on the cluster (4.9.0)
    //  that makes it impossible to use istringstream). The bug has been
    //  fixed in later versions of gcc.
    vector<stringstream*> content_no_burn_in;
    size_t line_count = 1;
    string line;
    while (!m_file.eof()) {
      getline(m_file, line);
      if (line_count > burn_in && line != "") {
        content_no_burn_in.push_back(new stringstream(line));
      }
      ++line_count;
    }
    m_file.close();

    // Check that after removing burnIn there are still enough iterations
    assert(content_no_burn_in.size() > n_params_sets);

    // Reshuffle lines
    gen.seed(seed);
    shuffle(content_no_burn_in.begin(), content_no_burn_in.end(), gen);

    // Reset tree output file
    ofstream treeFile(fTree);
    treeFile << "simID,time, infected,type_infected,infector,type_infector"
             << endl;
    treeFile.close();

    // Loop over params_set
    size_t simID;
    ofstream outputFile(out_file);
    for (size_t currP = 0; currP < n_params_sets; ++currP) {
      // Make sure we are reading from the beginning of stream
      content_no_burn_in[currP]->seekg(0, ios_base::beg);
      // Set parameters
      double to_skip, param_value;
      *(content_no_burn_in[currP]) >> to_skip; // LogLik
      for (size_t param_i = 0; param_i < nParams; ++param_i) {
        *(content_no_burn_in[currP]) >> param_value; // Parameters
        setParameter(param_i, param_value);
      }
      for (size_t currD = 0; currD < n_draws; ++currD) {
        simID = currP * n_draws + currD + 1;
        cout << "Simulation " << simID << endl;
        simulate(outputFile, simID, fromDate, toDate, fTree, ePalmi, eGalli);
      }
    }
    outputFile.close();

    // Clean up
    for (size_t curr = 0; curr < content_no_burn_in.size(); ++curr) {
      delete content_no_burn_in[curr];
    }
  }

  void getProbaInfectors(size_t cellJ, vector<size_t>& possibleInfectors,
      vector<double>& probaInfectors) {
    int infectionJ = cells[cellJ].dateInfection;
    size_t nInfectious = infectiousCellsAt[infectionJ].size();
    // nCells used as out of bound value (unknown infector)
    possibleInfectors = vector<size_t>(nInfectious + 1, nCells);
    probaInfectors = vector<double>(nInfectious + 1, 0.0);

    size_t curr = 0;
    double norm = 0.0;
    double curr_foi;
    for (size_t cellI : infectiousCellsAt[infectionJ]) {
      if (cellI == cellJ) continue;
      curr_foi = foi(cellJ, cellI, infectionJ, cells);
      possibleInfectors[curr] = cellI;
      probaInfectors[curr] = curr_foi;
      ++curr;
      norm += curr_foi;
    }
    // External infector
    curr_foi = foiExt(cellJ, infectionJ, cells);
    probaInfectors[nInfectious] = curr_foi;
    norm += curr_foi;

    // Normalize
    for (curr = 0; curr <= nInfectious; ++curr) {
      probaInfectors[curr] /= norm;
    }
  }

  void assignInfectors(vector<size_t>& infectees, vector<size_t>& infectors) {
    size_t curr = 0;
    for (size_t cellI : infectedCellsIDs) {
      infectees[curr] = cellI;
      vector<size_t> possibleInfectors;
      vector<double> probaInfectors;
      getProbaInfectors(cellI, possibleInfectors, probaInfectors);
      discrete_distribution<size_t> choice(probaInfectors.begin(),
        probaInfectors.end());
      infectors[curr] = possibleInfectors[choice(gen)];
      ++curr;
    }
  }

  void reconstructTransmissionTrees(string mcmc_file, string out_file,
      size_t burn_in, size_t n_params_sets, size_t n_draws,
      size_t seed = 0) {
    // Check that MCMC output file is there
    ifstream m_file(mcmc_file);
    if (!m_file.is_open()) {
      cout << "ERROR: file " << mcmc_file << " does not exist!" << endl;
      exit(1);
    }

    // Read in file (leaving out burn_in iterations)
    // NOTE: need to use stringstream pointers to run on cluster (there is
    //  a bug on the gcc version currently installed on the cluster (4.9.0)
    //  that makes it impossible to use istringstream). The bug has been
    //  fixed in later versions of gcc.
    vector<stringstream*> content_no_burn_in;
    size_t line_count = 1;
    string line;
    while (!m_file.eof()) {
      getline(m_file, line);
      if (line_count > burn_in && line != "") {
        content_no_burn_in.push_back(new stringstream(line));
      }
      ++line_count;
    }
    m_file.close();

    // Check that after removing burnIn there are still enough iterations
    assert(content_no_burn_in.size() > n_params_sets);

    // Reshuffle lines
    gen.seed(seed);
    shuffle(content_no_burn_in.begin(), content_no_burn_in.end(), gen);

    size_t nInfected = infectedCellsIDs.size();

    // Loop over params_set
    size_t simID;
    ofstream o_file(out_file);
    o_file << "Sim,Infectee,Infector,Distance,InfectionDate" << endl;
    for (size_t currP = 0; currP < n_params_sets; ++currP) {
      // Make sure we are reading from the beginning of stream
      content_no_burn_in[currP]->seekg(0, ios_base::beg);
      // Set parameters
      double to_skip, param_value;
      *(content_no_burn_in[currP]) >> to_skip; // LogLik
      for (size_t param_i = 0; param_i < nParams; ++param_i) {
        *(content_no_burn_in[currP]) >> param_value; // Parameters
        setParameter(param_i, param_value);
      }
      for (size_t currD = 0; currD < n_draws; ++currD) {
        simID = currP * n_draws + currD + 1;
        // nCells used as out of bound value
        vector<size_t> infectees(nInfected, nCells);
        vector<size_t> infectors(nInfected, nCells);
        assignInfectors(infectees, infectors);
        for (size_t curr = 0; curr < nInfected; ++curr) {
          o_file << simID << "," << cells[infectees[curr]].ID << ",";
          if (infectors[curr] == nCells) {
            o_file << "NA,NA,NA" << endl;
          } else {
            o_file << cells[infectors[curr]].ID << "," <<
              distanceMatrix[infectees[curr]][infectors[curr]] << ","
              << cells[infectees[curr]].dateInfection << endl;
          }
        }
      }
    }
    o_file.close();

    // Clean up
    for (size_t curr = 0; curr < content_no_burn_in.size(); ++curr) {
      delete content_no_burn_in[curr];
    }
  }

  void setParameter(size_t param_i, double value) {
    params[param_i] = value;
  }

  string toString() const {
    stringstream output;
    for (size_t curr = 0; curr < (nParams - 1); ++curr) {
      output << params[curr] << " ";
    }
    output << params[nParams - 1];
    return output.str();
  }
};


// -----------------------------------------------------------------------------
// Testing
// -----------------------------------------------------------------------------
void checkModel() {
  vector<int> switchPoints = {56, 76};
  SurvivalModel model("data.txt", 1, 1, 7, 7, 1, 5, 15.0, switchPoints);

  assert(std::abs(model.cutoffDistance - 15.0) < 1e-6);

  // Check that cellsWithinCutoffDistance is properly defined
  for (size_t celli = 0; celli < model.nCells; ++celli) {
    for (size_t cellj : model.cellsWithinCutoffDistance[celli]) {
      assert(model.distanceMatrix[celli][cellj] < model.cutoffDistance);
    }
  }

  // Check that loglik computation works as expected
  double ll1 = model.computeLogLik();
  double ll2 = model.computeLogLikV2();
  assert(std::abs(ll1 - ll2) < 1e-6);
}


void checkCode() {
  cout << "Checking model..." << endl;
  checkModel();
  cout << "All tests passed!" << endl;
}


// -----------------------------------------------------------------------------
// Helper functions to run MCMC
// -----------------------------------------------------------------------------
void runModel(string dataFile,
    int latentPeriod_Palmi, int latentPeriod_Galli,
    int incubationPeriod_Palmi, int incubationPeriod_Galli,
    int detectionDelay, int delayAsymptomatic,
    double cutoffDistance, const vector<int> &switchPoints, string oSuffix) {
  size_t seed = 131713;
  int toDate = 150;
  SurvivalModel model(dataFile, latentPeriod_Palmi, latentPeriod_Galli,
    incubationPeriod_Palmi, incubationPeriod_Galli, detectionDelay,
    delayAsymptomatic, cutoffDistance, switchPoints);

  MCMC mchain(&model, seed);
  mchain.outSuffix = oSuffix;

  // Set bounds
  for (size_t paramI = 0; paramI < model.nParams; ++paramI) {
    mchain.setBounds(paramI, 0.0, 1000.0);
  }

  // Find approximatly optimal random walk steps
  mchain.nIterations = 500;
  mchain.thinningFactor = 1;
  size_t maxReps = 50;
  double tolerance = 0.05;
  size_t repeat = 0;
  size_t decayStart = maxReps / 4;
  double delta0 = 2.0;
  double delta = delta0;
  while (true) {
    ++repeat;
    if (repeat >= decayStart) delta = delta0 / (repeat - decayStart + 1);
    cout << "Fine tuning " << repeat << endl;
    mchain.run();
    bool isOK = true;
    for (size_t i = 0; i < model.nParams; ++i) {
      if (mchain.rwSteps[i] > 0.0) {
        double curr_AR = mchain.acceptanceRates[i];
        double diff = curr_AR - 0.23;
        // Adjust random walk step
        mchain.rwSteps[i] *= (1.0 + delta * diff);
        isOK = isOK && (abs(diff) <= tolerance);
      }
    }
    if (isOK || repeat >= maxReps) break;
  }

  // Run with optimized random walk steps
  mchain.nIterations = 2500;
  mchain.thinningFactor = 10;
  mchain.run();

  // Tree reconstruction
  model.reconstructTransmissionTrees("output" + oSuffix + ".dat",
    "output" + oSuffix + "_trees.csv", 1000, 100, 5, seed);

  // Simulations
  model.simulateEpidemics("output" + oSuffix + ".dat",
    "output" + oSuffix + "_sims.txt", 1000, 100, 5, 1, toDate, seed);
}


void simulate_alternatives() {
  string dataFile = "data.txt";
  SurvivalModel model(dataFile, 1, 1, 7, 7, 1, 5, 15.0, {56, 76});

  int toDate = 150;
  string mcmc_file = "output_1_1_7_7_1_5_C15_S56_S76.dat";
  size_t seed = 956718921117;
  size_t burn_in = 1000;
  size_t n_params_sets = 100;
  size_t n_draws = 5;
  string sims_name;

  // Baseline
  sims_name = "sims_baseline.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);

  // No preventive culling
  model.preemptive_culling = false;
  sims_name = "sims_no_preventive_culling.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);

  // Varying preventive culling radius
  model.preemptive_culling = true;
  model.setCullingRadius0(1.0);
  model.setCullingRadius1(1.0);
  sims_name = "sims_preventive_culling_r_1.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);
  model.setCullingRadius0(2.0);
  model.setCullingRadius1(2.0);
  sims_name = "sims_preventive_culling_r_2.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);
  model.setCullingRadius0(3.0);
  model.setCullingRadius1(3.0);
  sims_name = "sims_preventive_culling_r_3.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);
  model.setCullingRadius0(4.0);
  model.setCullingRadius1(4.0);
  sims_name = "sims_preventive_culling_r_4.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);
  model.setCullingRadius0(5.0);
  model.setCullingRadius1(5.0);
  sims_name = "sims_preventive_culling_r_5.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);

  // Type of farms culled
  model.setCullingRadius0(1.0);
  model.setCullingRadius1(3.0);
  model.cullingTargets = 1;
  sims_name = "sims_preventive_culling_only_palmi.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);
  model.cullingTargets = 2;
  sims_name = "sims_preventive_culling_only_pag.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);

  // Culling delay for IPs
  model.cullingTargets = 0;
  model.rescaleGamma("culling", 2.0);
  sims_name = "sims_delay_culling_2.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);

  model.rescaleGamma("culling", 10.0);
  sims_name = "sims_delay_culling_10.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);

  // Preventive culling delay
  model.rescaleGamma("culling", 5.187740261066934);
  model.rescaleGamma("preventive_culling", 4.0);
  sims_name = "sims_delay_preventive_culling_4.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);
  model.rescaleGamma("preventive_culling", 16.0);
  sims_name = "sims_delay_preventive_culling_16.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);

  // Surveillance zones
  model.rescaleGamma("preventive_culling", 7.745806439112254);
  model.setSurveillanceRadius(5.0);
  sims_name = "sims_sz_5.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);
  model.setSurveillanceRadius(15.0);
  sims_name = "sims_sz_15.txt";
  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed);
}

// Simulate transmission tree
void simulate_trees(string fPrefix,
                    size_t n_params_sets = 100, size_t n_draws = 5,
                    string eFarm = "") {
  string dataFile = "data.txt";
  SurvivalModel model(dataFile, 1, 1, 7, 7, 1, 5, 15.0, {56, 76});

  int toDate = 150;
  string mcmc_file = "output_1_1_7_7_1_5_C15_S56_S76.dat";
  size_t seed = 956718921117;
  size_t burn_in = 1000;
  string sims_name = fPrefix + "_sims.txt";
  string trees_name = fPrefix + "_trees.txt";

  model.simulateEpidemics(mcmc_file, sims_name, burn_in, n_params_sets, n_draws,
    1, toDate, seed, trees_name, eFarm);
}

// -----------------------------------------------------------------------------
// Main
// -----------------------------------------------------------------------------
int main(int argc, char *argv[]) {
  // ---------------------------------------------------------------------------
  // To run few basic checks
  // ---------------------------------------------------------------------------
  // checkCode();

  // ---------------------------------------------------------------------------
  // To run the MCMC
  // ---------------------------------------------------------------------------
  /* necessary for the original model
  string dataFile = "data.txt";
  int latentPeriod_Palmi = 1;
  int latentPeriod_Galli = 1;
  int incubationPeriod_Palmi = 7;
  int incubationPeriod_Galli = 7;
  int detectionDelay = 1;
  int delayAsymptomatic = 5;
  double cutoffDistance = 15.0;
  vector<int> switchPoints = {56, 76};
  string mcmcSuff = "_1_1_7_7_1_5_C15_S56_S76";
  runModel(dataFile, latentPeriod_Palmi, latentPeriod_Galli,
    incubationPeriod_Palmi, incubationPeriod_Galli, detectionDelay,
    delayAsymptomatic, cutoffDistance, switchPoints, mcmcSuff);
  */

  // ---------------------------------------------------------------------------
  // To simulates the alternative control strategies described in the paper
  // ---------------------------------------------------------------------------
  // simulate_alternatives();
  string com_file = "", prefix = "baseline";
  int nparams = 100, ndraws = 5;
  if (argc != 4) {
    cerr << "usage: avian.exe baseline|com_file.txt nparams ndraws\n";
    exit(1);
  } else {
    com_file = string(argv[1]);
    nparams = atoi(argv[2]);
    ndraws = atoi(argv[3]);
  }
  if (com_file != "baseline") {
    int pos = com_file.find(".");
    prefix = com_file.substr(0, pos);
  } else {
    prefix = "baseline";
    com_file = "";
  }
  cout << "Running " << nparams * ndraws << " sims with output file "
       << prefix << "_sims.txt and " << prefix << "_trees.txt\n";
  simulate_trees(prefix, nparams, ndraws, com_file);

  return 0;
}
