#include "error.H"
#include "mixVoidFraction.H"
#include "addToRunTimeSelectionTable.H"
#include "mpi.h"

namespace Foam {

defineTypeNameAndDebug(mixVoidFraction, 0);

addToRunTimeSelectionTable(voidFractionModel, mixVoidFraction, dictionary);

/*! \brief 构造函数 */
mixVoidFraction::mixVoidFraction(const dictionary& dict,
                                 cfdemCloud& sm):
  voidFractionModel(dict, sm),
  propsDict_(dict.subDict(typeName + "Props")),
  verbose_(false),
  procBoundaryCorrection_(propsDict_.lookupOrDefault<Switch>("procBoundaryCorrection", false)),
  alphaMin_(readScalar(propsDict_.lookup("alphaMin"))),
  alphaLimited_(false),
  // scaleUpVol_(readScalar(propsDict_.lookup("scaleUpVol"))),
  tooMuch_(0.0),
  cfdemUseOnly_(false),
  hasInit_(false) {

  // 单个颗粒覆盖最多网格数量
  maxCellsNumPerFineParticle_ = propsDict_.lookupOrDefault<int>("maxCellsNumPerFineParticle", 29);
  maxCellsNumPerMiddleParticle_ = propsDict_.lookupOrDefault<int>("maxCellsNumPerMiddleParticle", 29);
  maxCellsNumPerCoarseParticle_ = propsDict_.lookupOrDefault<int>("maxCellsNumPerCoarseParticle", 1000);
  maxCellsPerParticle_ = propsDict_.lookupOrDefault<int>("maxCellsPerParticle", 1000);
  if (maxCellsPerParticle_ < maxCellsNumPerCoarseParticle_) {
    FatalError << "mixVoidFraction::mixVoidFraction(): maxCellsPerParticle_ must bigger than maxCellsNumPerCoarseParticle_"
      << abort(FatalError);
  }

  if (alphaMin_ > 1 || alphaMin_ < 0.01) {
    FatalError << "alphaMin shloud be > 1 and < 0.01." << abort(FatalError);
  }

  // 检查是否指定了 weight 和 porosity, 如果指定, 则读入 weight 和 porosity
  checkWeightNporosity(propsDict_);

  if (propsDict_.found("verbose")) { verbose_ = true; }

  if (propsDict_.found("cfdemUseOnly")) { cfdemUseOnly_ = readBool(propsDict_.lookup("cfdemUseOnly")); }

  // 检查 locateModel 是否符合 procBoundaryCorrection_ 的条件
  // need change when mixed!!!!!!!!!!!!!!!!!!
  // if (procBoundaryCorrection_) {
  //   if (!(particleCloud_.locateM().type() == "engineIB")) {
  //     FatalError << typeName << ": You are requesting procBoundaryCorrection, this requires the use of engineIB!\n"
  //       << abort(FatalError);
  //   }
  // } else {
  //   if (particleCloud_.locateM().type() == "engineIB") {
  //     FatalError << typeName << ": You are using engineIB, this requires using procBoundaryCorrection=true!\n"
  //       << abort(FatalError);
  //   }
  // }

  // if (scaleUpVol_ < 1) {
  //   FatalError << "scaleUpVol shloud be > 1." << abort(FatalError);
  // }

  // 计算标志点相对颗粒中心的坐标
  [this] (void) -> void {
    int idx = 0;
    // 设置中心标志点偏移
    for (int i = 0; i < 3; ++i) { offsets[idx][i] = 0.0; }
    idx += 1;
    // 计算两个半径, 这两个半径构成的球面将颗粒分为 1 : 14 : 14 这三个部分, 这三个部分分别是: 中心球体a, 球面层 b, 球面层 c
    // 第一个部分是一个中心与粒子中心重合的球体。该子颗粒的半径 r1 计算如下: V(r1) / V(R) = r1 ^ 3 / R ^ 3 = 1 / 29 ==> r1 = R * (1/29)^(1/3)
    // 其余体积是一个球面层, 必须分为两个等体积球面层。这两个球面层之间的径向边界位置很容易获得: V(r1) / V(R) = r1 ^ 3 / R ^ 3 = 15 / 29 ==> r1 = R * (15/29)^(1/3)
    double r1 = cbrt(1.0 / double(numberOfMarkerPoints));
    double r2 = cbrt(15.0 / double(numberOfMarkerPoints));

    // 将球面层 b 和 c 分别等体积分成 14 个部分
    // 第一球面层中各体积径向质心点的位置是: r(s,1) = 0.62761 * R
    // 第二球面层中各体积径向质心点的位置是: r(s,2) = 0.90853 * R
    scalar r[2] = { 0.75 * (r2*r2*r2*r2 - r1*r1*r1*r1) / (r2*r2*r2 - r1*r1*r1),
                    0.75 * (1.0 - r2*r2*r2*r2) / (1.0 - r2*r2*r2) };

    for (label ir = 0; ir < 2; ++ir) {
      // 通过球坐标系找到 8 个标志点
      for (scalar zeta = 0.25 * M_PI; zeta < 2.0 * M_PI; zeta += 0.5 * M_PI) {
        for(scalar theta = 0.25 * M_PI; theta < M_PI; theta += 0.5 * M_PI) {
          offsets[idx][0] = r[ir] * Foam::sin(theta) * Foam::cos(zeta);
          offsets[idx][1] = r[ir] * Foam::sin(theta) * Foam::sin(zeta);
          offsets[idx][2] = r[ir] * Foam::cos(theta);
          idx += 1;
        }
      }
      // 通过笛卡尔坐标系找到 6 个标志点(x y z 方向各有两个)
      for (int j = -1; j <= 1; j += 2) {
        offsets[idx][0] = r[ir] * static_cast<double>(j);
        offsets[idx][1] = 0.0;
        offsets[idx][2] = 0.0;
        idx += 1;

        offsets[idx][0] = 0.0;
        offsets[idx][1] = r[ir] * static_cast<double>(j);
        offsets[idx][2] = 0.0;
        idx += 1;

        offsets[idx][0] = 0.0;
        offsets[idx][1] = 0.0;
        offsets[idx][2] = r[ir] * static_cast<double>(j);
        idx += 1;
      }
    }
    return;
  } ();
}

/*! \brief 计算颗粒尺寸与其周围网格平均尺寸的比值, 并将颗粒索引按照颗粒尺寸归类 */
void mixVoidFraction::getDimensionRatios(std::vector<double>& dimensionRatios,
                                         std::vector<double>& globalDimensionRatios,
                                         std::vector<double>& sumCellsNumbers,
                                         std::vector<double>& sumCellsVolumes,
                                         const std::vector<double*>& centreCellIDs) {
  dimensionRatios.clear();
  globalDimensionRatios.clear();
  sumCellsNumbers.clear();
  sumCellsVolumes.clear();
  MPI_Barrier(MPI_COMM_WORLD);

  for (int index = 0; index < particleCloud_.numberOfParticles(); ++index) {
    // 获取颗粒中心所在网格编号
    label particleCenterCellID = *(centreCellIDs[index]);
    // 获取颗粒半径
    scalar radius = particleCloud_.radius(index);
    // 获取颗粒中心位置
    vector positionCenter = particleCloud_.position(index);

    if (particleCenterCellID >= 0) {  // particle found
      // 定义初始化模糊哈希集合
      labelHashSet initHashSett;
      // 构建初始化哈希集合
      buildLabelHashSetForDimensionRatio(index,
                                         positionCenter,
                                         particleCenterCellID,
                                         initHashSett, 1.0);
      // 计算颗粒周围网格的平均尺寸
      scalar Vmesh = 0.0;
      forAll (initHashSett, i) {
        label cellI = initHashSett.toc()[i];
        Vmesh += particleCloud_.mesh().V()[cellI];
      }
      // 计算颗粒尺寸与颗粒中心所在网格尺寸的比值
      double ratio = pow(Vmesh / static_cast<double>(initHashSett.size()), 1.0 / 3.0) / (2.0 * radius);
      dimensionRatios.push_back(ratio);
      sumCellsNumbers.push_back(static_cast<double>(initHashSett.size()));
      sumCellsVolumes.push_back(static_cast<double>(Vmesh));
    } else {
      // 如果在当前处理器没有搜索到颗粒覆盖的网格, 则置为 -1.0
      dimensionRatios.push_back(-1);
      sumCellsNumbers.push_back(-1);
      sumCellsVolumes.push_back(-1);
    }
    // 对所有颗粒都置为 -1
    globalDimensionRatios.push_back(-1);
  }  // End of index

  // 这里必须调用 particleCloud_.mesh().C(), 否则在下面的 MPI_Barrier 中进程会阻塞
  Pout << "Mesh number in current Proc: " << particleCloud_.mesh().C().size() << endl;
  MPI_Barrier(MPI_COMM_WORLD);

  // 主节点汇总其他节点的 dimension ratio 等信息
  int numberOfParticles = particleCloud_.numberOfParticles();
  int numProc, myProc;
  MPI_Comm_size(MPI_COMM_WORLD, &numProc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myProc);
  int globalDimensionRatios_tag = 100;
  int sumCellsNumbers_tag = 102;
  int sumCellsVolumes_tag = 103;

  if (myProc != 0) {
    // 不是主节点
    MPI_Request myRequest[3];
    MPI_Status myStatus[3];
    // 发送 sumCellsNumbers 给主节点( MPI_Isend 非阻塞)
    MPI_Isend(sumCellsNumbers.data(), numberOfParticles, MPI_DOUBLE, 0,
              sumCellsNumbers_tag, MPI_COMM_WORLD, &myRequest[0]);
    MPI_Isend(sumCellsVolumes.data(), numberOfParticles, MPI_DOUBLE, 0,
              sumCellsVolumes_tag, MPI_COMM_WORLD, &myRequest[1]);
    // 从主节点接收 globalDimensionRatios
    MPI_Irecv(globalDimensionRatios.data(), numberOfParticles, MPI_DOUBLE, 0,
              globalDimensionRatios_tag, MPI_COMM_WORLD, &myRequest[2]);
    MPI_Wait(myRequest + 2, myStatus + 2);
  }

  if (myProc == 0) {
    std::vector<MPI_Request> myRequest_vec(numProc);
    std::vector<MPI_Status> myStatus_vec(numProc);
    std::vector<double> sumCellsNumbers_vec(numProc * numberOfParticles, -1);
    std::vector<double> sumCellsVolumes_vec(numProc * numberOfParticles, -1);

    // 如果是主节点, 则接收其他节点的 sumCellsNumbers 信息
    for (int inode = 1; inode < numProc; ++inode) {
      MPI_Irecv(sumCellsNumbers_vec.data() + inode * numberOfParticles,
                numberOfParticles, MPI_DOUBLE, inode, sumCellsNumbers_tag,
                MPI_COMM_WORLD, myRequest_vec.data() + inode);
    }  // End of loop sub node
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, myRequest_vec.data() + 1, myStatus_vec.data() + 1);

    // 如果是主节点, 则接收其他节点的 sumCellsVolumes 信息
    for (int inode = 1; inode < numProc; ++inode) {
      MPI_Irecv(sumCellsVolumes_vec.data() + inode * numberOfParticles,
                numberOfParticles, MPI_DOUBLE, inode, sumCellsVolumes_tag,
                MPI_COMM_WORLD, myRequest_vec.data() + inode);
    }  // End of loop sub node
    // 主节点等待 Irecv 执行完成
    MPI_Waitall(numProc - 1, myRequest_vec.data() + 1, myStatus_vec.data() + 1);

    // 由主节点计算 globalDimensionRatios
    for (int i = 0; i < numberOfParticles; ++i) {
      double Nmesh = 0, Vmesh = 0;
      double radius = particleCloud_.radius(i);
      for (int j = 0; j < numProc; ++j) {
        if (j == 0) {
          Nmesh += std::max(sumCellsNumbers[i], 0.);
          Vmesh += std::max(sumCellsVolumes[i], 0.);
        }
        Nmesh += std::max(sumCellsNumbers_vec[i + j * numberOfParticles], 0.);
        Vmesh += std::max(sumCellsVolumes_vec[i + j * numberOfParticles], 0.);
      }
      double globalRatio = pow(Vmesh / Nmesh, 1.0 / 3.0) / (2.0 * radius);
      globalDimensionRatios[i] = globalRatio;
    }

    for (int inode = 1; inode < numProc; ++inode) {
      MPI_Isend(globalDimensionRatios.data(), numberOfParticles, MPI_DOUBLE, inode,
                globalDimensionRatios_tag, MPI_COMM_WORLD, myRequest_vec.data() + inode);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

/*!
 * \brief 计算被 middle 颗粒影响到的网格编号
 * \note cellAffectIDs 用于计算 middle particle 受到的流体阻力, 搜索范围默认是 4 倍颗粒半径
 */
void mixVoidFraction::getAffectedCellIDs(const std::vector<int>& middleParticleIndexs,
                                         std::vector<std::vector<int> >& affectedCellIDs,
                                         const scalar& searchScale) const {
  if (searchScale < 2 && searchScale > 20) {
    FatalError << "mixVoidFraction::getAffectedCellIDs(): searchScale " << searchScale << " not in range "
      << "[2, 100]" << abort(FatalError);
  }

  // 初始化 affectedCellIDs
  affectedCellIDs.clear();

  for (const int& index : middleParticleIndexs) {
    // 获取颗粒中心所在网格编号
    label particleCenterCellID = particleCloud_.cellIDs()[index][0];
    // 获取颗粒中心位置
    vector positionCenter = particleCloud_.position(index);

    if (particleCenterCellID >= 0) {  // particle found
      // 定义被颗粒影响的网格的哈希集合
      labelHashSet affectedHashSett;
      buildLabelHashSetForDimensionRatio(index, positionCenter,
                                         particleCenterCellID,
                                         affectedHashSett,
                                         searchScale);
      // 保存搜索结果
      affectedCellIDs.push_back(std::vector<int>());
      auto iter = affectedCellIDs.end();
      for (int i = 0; i < affectedHashSett.size(); ++i) {
        iter->push_back(affectedHashSett.toc()[i]);
      }
    } else {
      affectedCellIDs.push_back(std::vector<int>());
    }  // End of particle found
  }  // End of loop middleParticleIndexs
}

/*!
 * \brief (模糊搜索!!)构建颗粒覆盖的所有网格的哈希集合
 * \note 设置为递归函数，通过哈希器将网格编号转换为哈希值，并存入 set 中以便于搜索
 * \param index     <[in] 颗粒索引
 * \param position  <[in] 颗粒中心位置
 * \param cellID    <[in] 递归循环中要检索网格编号
 * \param hashSett  <[in, out] 需要构建的哈希集
 * \param scale     <[in] 颗粒半径扩大系数
 */
void mixVoidFraction::buildLabelHashSetForDimensionRatio(const int index,
                                                         const vector& position,
                                                         const label& cellID,
                                                         labelHashSet& hashSett,
                                                         const double& scale) const {
  if (cellID < 0) { return; }

  hashSett.insert(cellID);
  // 获取 cellID 网格的所有 neighbour cell 的链表
  const labelList& nc = particleCloud_.mesh().cellCells()[cellID];

  // 遍历链表
  forAll(nc, i) {
    // 获取相邻网格索引, 以及网格中心坐标
    double neighbour = nc[i];

    if (neighbour >= 0) {  // 判断是否搜索到边界
      vector cellCentrePosition = particleCloud_.mesh().C()[neighbour];

      if (false == hashSett.found(neighbour)) {
        if (pointInParticle(index, position, cellCentrePosition, scale) < 0.0 ) {  // 如果 neighbour 网格中心在颗粒中, 并且在哈希集合中没有 neighbour 网格的索引
          // 以 neighbour 为中心递归构建哈希集合
          buildLabelHashSetForDimensionRatio(index, position, neighbour, hashSett, scale);
        }
      }
    }
  }
}

/*!
 * \brief 设置空隙率
 * \note for other solvers used
 */
void mixVoidFraction::setvoidFraction(double** const& mask,
                                      double**& voidfractions,
                                      double**& particleWeights,
                                      double**& particleVolumes,
                                      double**& particleV) const {
  if (cfdemUseOnly_) {
    reAllocArrays(particleCloud_.numberOfParticles());
  } else {
    reAllocArrays();
  }

  if (particleCloud_.usedForSolverIB()) {
    voidfractionNext_ == dimensionedScalar("one", voidfractionNext_.dimensions(), 1.);
  }

  for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
    // skip this particle if not correct type
    if (!checkParticleType(index)) { continue; }

    // 对索引为 index 的颗粒初始化
    for (int subcell = 0; subcell < cellsPerParticle_[index][0]; subcell++) {
      particleWeights[index][subcell] = 0.;
      particleVolumes[index][subcell] = 0.;
    }
    particleV[index][0] = 0.;
    cellsPerParticle_[index][0] = 1.;

    if (particleCloud_.usedForSolverIB()) {
      // 设置索引为 index 的颗粒的流体体积分数
      setVolumeFractionForSingleParticle(index);
    } else if (particleCloud_.usedForSolverPiso()) {
      // 设置索引为 index 的颗粒的空隙率
      setvoidFractionForSingleParticle(index,
                                       particleWeights,
                                       particleVolumes,
                                       particleV);
    }
  }  // End of loop all particles

  if (particleCloud_.usedForSolverIB()) {
    // bring voidfraction from Eulerian Field to particle array
    for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
      // 遍历颗粒覆盖的所有网格
      for (int subcell = 0; subcell < cellsPerParticle_[index][0]; subcell++) {
        // 获取网格编号
        label cellID = particleCloud_.cellIDs()[index][subcell];

        if (cellID >= 0) {
          voidfractions[index][subcell] = voidfractionNext_[cellID];
        } else {
          voidfractions[index][subcell] = -1.0;
        }
      }  // End of subcell
    }
  } else if (particleCloud_.usedForSolverPiso()) {
    voidfractionNext_.correctBoundaryConditions();
    // reset counter of lost volume
    if (verbose_) { Pout << "Total particle volume neglected: " << tooMuch_<< endl; }
    tooMuch_ = 0.;

    // bring voidfraction from Eulerian Field to particle array
    for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
      // 遍历颗粒覆盖的所有网格
      for (int subcell = 0; subcell < cellsPerParticle_[index][0]; subcell++) {
        // 获取网格编号
        label cellID = particleCloud_.cellIDs()[index][subcell];
        if (cellID >= 0) {
          voidfractions[index][subcell] = voidfractionNext_[cellID];
        } else {
          voidfractions[index][subcell] = -1.;
        }
      }  // End of subcell
    }
  }  // End of particleCloud_.usedForSolverIB()
}

/*! \brief 设置空隙率以及颗粒体积分数前必须调用该初始化函数 */
void mixVoidFraction::voidFractionModelInit(double**& particleWeights,
                                            double**& particleVolumes,
                                            double**& particleV,
                                            const std::vector<double>& dimensionRatios) const {
  // 如果已经初始化过则直接返回
  if (hasInit_) { return; }

  // (1) 重新分配 cellsPerParticle_
  if (cfdemUseOnly_) {
    reAllocArrays(particleCloud_.numberOfParticles());
  } else {
    reAllocArrays();
  }

  // (2) 重置 tooMuch_
  tooMuch_ = 0;

  for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
    // (3) 所有颗粒将 cellsPerParticle_ 重置为 1
    cellsPerParticle_[index][0] = 1.;

    double ratio = dimensionRatios[index];

    // (4) fine and middle 颗粒重置
    if (particleCloud_.checkFAndMParticle(ratio)) {
      for (int subcell = 0; subcell < std::max(maxCellsNumPerFineParticle_, maxCellsNumPerMiddleParticle_); subcell++) {
        particleWeights[index][subcell] = 0.;
        particleVolumes[index][subcell] = 0.;
      }
      particleV[index][0] = 0.;
    }  // End of fine or middle particles
  }  // End of index

  // (5) 如果在 evolve 函数中没有调用 resetVoidFractions() or resetVolumeFractions()，那么这里需要调用
  hasInit_ = true;
}

/*!
 * \brief 设置空隙率以及流体体积分数
 * \note for cfdemSolverMix used
 */
void mixVoidFraction::setVoidFractionAndVolumeFraction(double** const& mask,
                                                       double**& volumefractions,
                                                       double**& voidfractions,
                                                       double**& particleWeights,
                                                       double**& particleVolumes,
                                                       double**& particleV,
                                                       const std::vector<double>& dimensionRatios) const {
  if (!hasInit_) {
    FatalError << "mixVoidFraction::setvoidFraction(): hasInit_ == false, please invoke mixVoidFraction() function"
      << abort(FatalError);
  }

  for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
    // skip this particle if not correct type
    if (!checkParticleType(index)) { continue; }

    double ratio = dimensionRatios[index];
    if (particleCloud_.needSetFieldForCoarseParticle(index, false, dimensionRatios)) {
      setVolumeFractionForSingleParticle(index);
    } else {
      // 如果颗粒是 fine or middle，设置索引为 index 的颗粒的空隙率
      setvoidFractionForSingleParticle(index,
                                       particleWeights,
                                       particleVolumes,
                                       particleV);
    }
  }  // End of loop all particles

  voidfractionNext_.correctBoundaryConditions();
#if 0
  // 如果使用 dynamicRefineMesh, 则在此处不能修正边界条件
  volumefractionNext_.correctBoundaryConditions();
#endif

  if (verbose_) { Pout << "Total particle volume neglected: " << tooMuch_<< endl; }

  // 将空隙率场和流体体积分数场从欧拉场转换到拉格朗日场
  int fmWidth = std::max(maxCellsNumPerFineParticle_, maxCellsNumPerMiddleParticle_);
  for (int index = 0; index < particleCloud_.numberOfParticles(); index++) {
    double ratio = dimensionRatios[index];
    if (particleCloud_.needSetFieldForCoarseParticle(index, false, dimensionRatios)) {
      // 将 coarse particle 的 voidfraction 全部设置为 -1.0
      for (int subcell = 0; subcell < fmWidth; subcell++) {
        voidfractions[index][subcell] = -1.0;
      }
      // 遍历颗粒覆盖的所有网格
      for (int subcell = 0; subcell < maxCellsNumPerCoarseParticle_; subcell++) {
        if (subcell < cellsPerParticle_[index][0]) {
          // 获取网格编号
          label cellID = particleCloud_.cellIDs()[index][subcell];
          if (cellID >= 0) {
            volumefractions[index][subcell] = volumefractionNext_[cellID];
          } else {
            volumefractions[index][subcell] = -1.0;
          }
        } else {
          volumefractions[index][subcell] = -1.0;
        }
      }
    } else {
      // 遍历颗粒覆盖的所有网格
      for (int subcell = 0; subcell < fmWidth; subcell++) {
        // 将 fine and middle particle 的 volumefraction 全部设置为 -1.0
        volumefractions[index][subcell] = -1.0;
        if (subcell < cellsPerParticle_[index][0]) {
          // 获取网格编号
          label cellID = particleCloud_.cellIDs()[index][subcell];
          if (cellID >= 0) {
            voidfractions[index][subcell] = voidfractionNext_[cellID];
          } else {
            voidfractions[index][subcell] = -1.0;
          }
        } else {
          voidfractions[index][subcell] = -1.0;
        }
      }
    }
  } // End of loop all particles

  // Note: change hasInit_ to false
  hasInit_ = false;
}

/*! \brief 设置索引为 index 的单个颗粒的空隙率 */
void mixVoidFraction::setvoidFractionForSingleParticle(const int index,
                                                       double**& particleWeights,
                                                       double**& particleVolumes,
                                                       double**& particleV) const {
  if (index < 0 || index >= particleCloud_.numberOfParticles()) {
    FatalError << "mixVoidFraction.C::setvoidFractionForSingleParticle(): " << "index " << index
      << "out of boundary [0, " << particleCloud_.numberOfParticles() - 1 << "]" << abort(FatalError);
  }

  scalar cellVol(0.);
  scalar scaleVol = weight();
  scalar scaleRadius = cbrt(porosity());
  // const boundBox& globalBb = particleCloud_.mesh().bounds();

  // 获取颗粒 index 中心坐标, 颗粒中心所在网格编号, 以及颗粒半径
  vector position = particleCloud_.position(index);
  label cellID = particleCloud_.cellIDs()[index][0];
  scalar radius = particleCloud_.radius(index);
  // 计算颗粒体积 = 颗粒体积 * 颗粒体积因数
  scalar volume = 4.188790205 * radius * radius * radius * scaleVol;
  // 先计算体积, 然后再用 scaleRadius 乘以半径, 以保证在计算空隙率时候颗粒尺寸不变
  radius *= scaleRadius;

  // 定义 variables for sub-search
  int cellsSet = 0;

  // label cellWithCenter(-1);
  // if (procBoundaryCorrection_) {
  //   // 如果修正处理器边界, 则重新定位颗粒中心所在的网格坐标
  //   cellWithCenter = particleCloud_.locateM().findSingleCell(position, cellID);
  //   particleCloud_.cellIDs()[index][0] = cellWithCenter;
  // }

  if (cellID >= 0) {  // particle centre is in domain
    cellVol = particleCloud_.mesh().V()[cellID];

    // 遍历所有标志点(除了颗粒中心处的标志点)
    for (int i = 1; i < numberOfMarkerPoints; i++) {

      // 获取标志点的绝对位置
      vector subPosition = position + radius * offsets[i];

      // // 检查周期性边界
      // if (particleCloud_.checkPeriodicCells()) {
      //   for (int iDir = 0; iDir < 3; iDir++) {
      //     // 如果标志点在计算区域外部, 那么获取在计算区域内部的标志点
      //     if (subPosition[iDir] > globalBb.max()[iDir]) {
      //       subPosition[iDir] -= globalBb.max()[iDir] - globalBb.min()[iDir];
      //     } else if (subPosition[iDir] < globalBb.min()[iDir]) {
      //       subPosition[iDir] += globalBb.max()[iDir] - globalBb.min()[iDir];
      //     }
      //   }
      // }

      // 根据修正后的标志点坐标定位标志点所在网格的索引
      label partCellId = particleCloud_.locateM().findSingleCell(subPosition, cellID);

      if (partCellId >= 0) {  // 如果标志点在 domain 中, 则更新空隙率
          // 获取 partCellId 网格体积
          scalar partCellVol = particleCloud_.mesh().V()[partCellId];
          // 获取颗粒体积的 1 / 29
          scalar particleVolume = volume / static_cast<scalar>(numberOfMarkerPoints);
          // 计算空隙率
          scalar newAlpha = voidfractionNext_[partCellId] - particleVolume / partCellVol;

          // 更新空隙率场
          if (newAlpha > alphaMin_) {
            voidfractionNext_[partCellId] = newAlpha;
          } else {
            // 如果空隙率低于最低空隙率, 则直接赋值为 alphaMin_, 并累加损失的颗粒体积
            voidfractionNext_[partCellId] = alphaMin_;
            tooMuch_ += (alphaMin_ - newAlpha) * partCellVol;
          }

          // 统计标志点数目
          cellsSet += 1;

          // 判断是否需要更新颗粒覆盖的网格索引
          bool createNew = true;
          label storeInIndex = 0;
          for (int j = 0; j < cellsPerParticle_[index][0]; ++j) {
            if (partCellId == particleCloud_.cellIDs()[index][j]) {
              // 如果 partCellId 已经在 cellIDs[index] 中, 则不需要在 cellIDs 中创建新的索引
              storeInIndex = j;
              createNew = false;
              break;
            }
          }
          if (storeInIndex >= maxCellsNumPerFineParticle_ || storeInIndex >= maxCellsNumPerMiddleParticle_) {
            FatalError << "mixVoidFraction::setvoidFractionForSingleParticle(): storeInIndex is out of range"
              << abort(FatalError);
          }
          // 如果需要创建新的索引
          if (createNew) {
            cellsPerParticle_[index][0] ++;
            storeInIndex = cellsPerParticle_[index][0] - 1;
            particleCloud_.cellIDs()[index][storeInIndex] = partCellId;
          }
          particleWeights[index][storeInIndex] += 1.0 / static_cast<scalar>(numberOfMarkerPoints);
          particleVolumes[index][storeInIndex] += particleVolume;
          particleV[index][0] += particleVolume;
      }  // End of partCellId >= 0
    }  // End of marker point loop

    if (cellsSet + 1 > numberOfMarkerPoints) {
      Info << "mixVoidFraction::setvoidFractionForSingleParticle(): wrong cellsSet = " << cellsSet << endl;
    }

    // 处理颗粒中心处的标志点, 将丢失的标志点对空隙率的影响作用在颗粒中心处的标志点上
    // if (!procBoundaryCorrection_) {
    if (1) {
      // set source for particle center
      scalar centreWeight = (numberOfMarkerPoints - cellsSet) * (1. / numberOfMarkerPoints);
      // update voidfraction for each particle read
      scalar newAlpha = voidfractionNext_[cellID] - volume * centreWeight / cellVol;
      if (newAlpha > alphaMin_) {
        voidfractionNext_[cellID] = newAlpha;
      } else {
        voidfractionNext_[cellID] = alphaMin_;
        tooMuch_ += (alphaMin_ - newAlpha) * cellVol;
      }
      // 设置颗粒中心标志点的影响系数
      particleWeights[index][0] += centreWeight;
      // store particleVolume for each particle
      particleVolumes[index][0] += volume * centreWeight;
      particleV[index][0] += volume * centreWeight;
    }
  }  // End of cellID >= 0
}

/*! \brief 设置索引为 index 的单个颗粒的体积分数场 */
void mixVoidFraction::setVolumeFractionForSingleParticle(const int index) const {

  if (index < 0 || index >= particleCloud_.numberOfParticles()) {
    FatalError << "mixVoidFraction.C::setVolumeFractionForSingleParticle(): " << "index " << index
      << "out of boundary [0, " << particleCloud_.numberOfParticles() - 1 << "]" << abort(FatalError);
  }
  // 设置需要计算的空隙率场
  volScalarField& fraction = particleCloud_.usedForSolverIB() ? voidfractionNext_ : volumefractionNext_;

  // 获取颗粒 index 中心坐标, 颗粒中心所在网格编号, 以及颗粒半径
  vector positionCenter = particleCloud_.position(index);
  label particleCenterCellID = particleCloud_.cellIDs()[index][0];
  // const boundBox& globalBb = particleCloud_.mesh().bounds();

  if (particleCenterCellID >= 0) {  // particle centre is in domain

    // 判断网格中心是否在颗粒中
    vector cellCentrePosition = particleCloud_.mesh().C()[particleCenterCellID];
    scalar fc = pointInParticle(index, positionCenter, cellCentrePosition);

    // 定义距离网格中心 cellCentrePosition 距离最小的颗粒坐标(可能是镜像颗粒的坐标, 也可能是原颗粒的坐标)
    // vector minPeriodicParticlePos = positionCenter;
    // if (particleCloud_.checkPeriodicCells()) {
    //   fc = minPeriodicDistance(index,                   // <[in] 颗粒的编号
    //                            cellCentrePosition,      // <[in] 网格中心的位置坐标
    //                            positionCenter,          // <[in] 原颗粒的坐标
    //                            globalBb,                // <[in] Boundary Box
    //                            minPeriodicParticlePos,  // <<[in, out] 距离最小颗粒坐标
    //                            particleCloud_.wall_periodicityCheckRange());
    // }
    // 计算网格 particleCenterCellID 的等效半径
    scalar corona = 0.5 * sqrt(3.0) * cbrt(particleCloud_.mesh().V()[particleCenterCellID]);

    // 获取网格 particleCenterCellID 的 corona point
    vector coronaPoint = getCoronaPointPosition(positionCenter, cellCentrePosition, corona);

    if (pointInParticle(index, positionCenter, coronaPoint) < 0.0) {
      // 如果 coronaPoint 在颗粒中, 则认为整个网格在颗粒中
      if (particleCloud_.usedForSolverIB()) {
        // 当求解器为 cfdemSolverIB 时, 在其他模型中使用的是 voidfractionNext_, 所以直接赋值给 voidfractionNext_
        voidfractionNext_[particleCenterCellID] = 0.0;
      } else {
        volumefractionNext_[particleCenterCellID] = 0.0;
      }
    } else {
      // 如果 coronaPoint 不在颗粒中, 则需要遍历网格的所有角点, 判断角点与网格中心是否在颗粒中
      const labelList& vertexPoints = particleCloud_.mesh().cellPoints()[particleCenterCellID];
      double ratio = 0.125;
      forAll (vertexPoints, i) {
        // 获取第 i 角点坐标
        vector vertexPosition = particleCloud_.mesh().points()[vertexPoints[i]];

        // 判断角点是否在颗粒中
        scalar fv = pointInParticle(index, positionCenter, vertexPosition);

        if (fc < 0.0 && fv < 0.0) {
          // 网格中心在颗粒中, 角点也在颗粒中
          fraction[particleCenterCellID] = 0.0;
        } else if (fc < 0.0 && fv > 0.0) {
          // 网格中心在颗粒中, 角点不在颗粒中
          // 计算角点 vertexPosition 对体积分数的影响系数 lambda
          scalar lambda = segmentParticleIntersection(index,
                                                      positionCenter,
                                                      cellCentrePosition,
                                                      vertexPosition);
          fraction[particleCenterCellID] -= ratio * lambda;
        } else if (fc > 0.0 && fv < 0.0) {
          // 网格中心不在颗粒中, 角点在颗粒中
          // 计算角点 vertexPosition 对体积分数的影响系数 lambda
          scalar lambda = segmentParticleIntersection(index,
                                                      positionCenter,
                                                      vertexPosition,
                                                      cellCentrePosition);
          if (particleCloud_.usedForSolverIB()) {
            voidfractionNext_[particleCenterCellID] -= ratio * lambda;
          } else {
            volumefractionNext_[particleCenterCellID] -= ratio * lambda;
          }
        }
      }  // End of loop of vertexPoints
    }
    // 颗粒中心所在网格的体积分数已经计算完成, 下面开始递归构建相邻网格
    labelHashSet hashSett;
    // 构建哈希集合
    buildLabelHashSetForVolumeFractions(index,
                                        positionCenter,
                                        particleCenterCellID,
                                        true,
                                        hashSett,
                                        fraction);
    scalar hashSetLength = hashSett.size();
    if (hashSetLength > maxCellsNumPerCoarseParticle_) {  // 如果集合中元素个数大于颗粒覆盖网格数限制
      FatalError << "Big particle found " << hashSetLength << " cells more than permittd maximun number of cells per paticle "
        << maxCellsNumPerCoarseParticle_ << abort(FatalError);
    } else if (hashSetLength > 0) {

      // 将颗粒覆盖的网格数保存到 cellsPerParticle_[index][0]
      cellsPerParticle_[index][0] = hashSetLength;

      // 去除已经存入 cellIDs()[index][0] 处的 particleCenterCellID
      hashSett.erase(particleCenterCellID);

      // 保存颗粒覆盖的所有网格编号
      for (label i = 0; i < hashSetLength - 1; ++i) {
        particleCloud_.cellIDs()[index][i + 1] = hashSett.toc()[i];
      }
    }
  }  // End of particleCenterCellID >= 0
}

/*!
 * \brief 构建颗粒覆盖的所有网格的哈希集合
 * \note 设置为递归函数,  通过哈希器将网格编号转换为哈希值,  并存入 set 中以便于搜索
 * \param index          <[in] 颗粒索引
 * \param position       <[in] 颗粒中心位置
 * \param cellID         <[in] 递归循环中要检索网格编号
 * \param initialInsert  <[in] 是否将 cellID 插入 hashSett
 * \param hashSett       <[in, out] 需要构建的哈希集
 * \param fraction       <[in, out] 需要设置的空隙率场
 */
void mixVoidFraction::buildLabelHashSetForVolumeFractions(const int& index,
                                                          const vector& position,
                                                          const label& cellID,
                                                          const bool& initialInsert,
                                                          labelHashSet& hashSett,
                                                          volScalarField& fraction) const {
  if (initialInsert) {
    hashSett.insert(cellID);
  }

  // 获取 cellID 网格的所有 neighbour cell 的链表
  const labelList& nc = particleCloud_.mesh().cellCells()[cellID];
  forAll (nc, i) {
    // 获取相邻网格索引,  以及网格中心坐标
    label neighbour = nc[i];
    if (neighbour < 0) {
      continue;
    } else {
      // 获取相邻网格中心坐标
      vector cellCentrePosition = particleCloud_.mesh().C()[neighbour];

      // 判断相邻网格中心是否在颗粒中
      scalar fc = pointInParticle(index, position, cellCentrePosition);

      // 计算相邻网格的等效半径
      scalar coronaRaidus = 0.5 * sqrt(3.0) * cbrt(particleCloud_.mesh().V()[neighbour]);

      // 获取 corona point
      vector coronaPoint = getCoronaPointPosition(position, cellCentrePosition, coronaRaidus);

      if (!hashSett.found(neighbour)) {
        // 如果在哈希集合中没有插入 neighbour 网格, 则需要计算 neighbour 网格的体积分数
        if (pointInParticle(index, position, coronaPoint) < 0.0) {
          // 如果相邻网格的 coronaPoint 在颗粒中, 则说明该网格完全被颗粒覆盖
          fraction[neighbour] = 0.0;

          // 以相邻网格为中心继续递归构建哈希集合
          buildLabelHashSetForVolumeFractions(index, position, neighbour,
                                              true, hashSett, fraction);
        } else {
          // 如果相邻网格的 coronaPoint 不在颗粒中, 则需要遍历该网格的所有角点
          // 定义单个角点对空隙率的影响率
          double ratio = 0.125;

          // 体积分数初始值
          scalar scale = 1.0;

          // 获取 neighbour 网格的角点集合
          const labelList& vertexPoints = particleCloud_.mesh().cellPoints()[neighbour];

          forAll (vertexPoints, j) {
            // 获取角点坐标
            vector vertexPosition = particleCloud_.mesh().points()[vertexPoints[j]];

            // 判断角点是否在颗粒中
            scalar fv = pointInParticle(index, position, vertexPosition);

            if (fc < 0.0 && fv < 0.0) {  // 如果网格 neighbour 中心在颗粒中, 角点 j 也在颗粒中
              // 体积分数直接减去 0.125
              scale -= ratio;
            } else if (fc < 0.0 && fv >= 0.0) {  // 如果网格 neighbour 中心在颗粒中, 角点 j 不在颗粒中
              // 计算角点 vertexPoints[j] 对空隙率的影响系数 lambda
              scalar lambda = segmentParticleIntersection(index,
                                                          position,
                                                          cellCentrePosition,
                                                          vertexPosition);
              // 体积分数减去 ratio * lambda
              scale -= ratio * lambda;
            } else if (fc >= 0.0 && fv < 0.0) {  // 如果网格 neighbour 中心不在颗粒中, 角点 j 在颗粒中
              // 计算角点 vertexPoints[j] 对空隙率的影响系数 lambda
              scalar lambda = segmentParticleIntersection(index,
                                                          position,
                                                          vertexPosition,
                                                          cellCentrePosition);
              // 体积分数减去 ratio * lambda
              scale -= ratio * lambda;
            }
          }  // End of loop vertexPoints

          // 保证体积分数 >= 0
          scale = scale >= 0.0 ? scale : 0.0;

          if (fabs(fraction[neighbour] - 1.0) < 1e-15) {
            // 如果 neighbour 网格的体积分数为 1.0, 则说明第一次遍历到该网格, 可以直接赋值
            fraction[neighbour] = scale;
          } else {
            // 如果 neighbour 网格的体积分数不为 1.0, 则说明在计算其他颗粒时候, 已经遍历到该网格
            fraction[neighbour] -= (1.0 - scale);
            // 保证体积分数 >= 0
            fraction[neighbour] = fraction[neighbour] >= 0.0 ? fraction[neighbour] : 0.0;
          }
          if (!(fabs(scale - 1.0) < 1e-15)) {
            // 如果体积分数不为 1.0, 则说明该 neighbour 需要递归循环构建哈希集合
            buildLabelHashSetForVolumeFractions(index, position, neighbour,
                                                true, hashSett, fraction);
          }
        }  // End of corona point not in particle
      }  // End of hashSett.found(neighbour)
    }  // End of neighbour >= 0
  }  // End of loop neighbour cell
}

/*!
 * \brief 获取 Corona Point
 * \note Corona Point 是在以网格中心为中心, 半径为 corona 的球面上, 距离颗粒中心最远的点
 *       其中, 半径 corona = 0.5 * sqrt(3) * 网格等效半径
 *       网格等效半径 = pow(cellVolume, 1.0 / 3.0)
 * \param position            <[in] 指定颗粒中心
 * \param cellCentrePosition  <[in] 指定网格中心
 * \param corona              <[in] 指定网格的等效半径
 */
vector mixVoidFraction::getCoronaPointPosition(vector const& position,
                                               vector const& cellCentrePosition,
                                               scalar const corona) const {
  // 计算网格中心到颗粒中心的距离
  scalar centreDist = mag(cellCentrePosition - position);
  vector coronaPoint = cellCentrePosition;
  if(centreDist > 0.0){
    coronaPoint = cellCentrePosition + (cellCentrePosition - position) * (corona / centreDist); 
    return coronaPoint; 
  }
  return coronaPoint;
}

/*!
 * \brief 计算距离系数
 * \note 对任意一个网格, 如果网格中心 c 在颗粒内部, 但是它的某个角点 p 不在颗粒内部,
 *       则计算 c 与 p 的连线与颗粒表面的交点 i 到网格中心 c 的距离, 即求解 x 的二元一次方程
 *       (x * (vector_p - vector_c) - vector_particle) & 
 *       (x * (vector_p - vector_c) - vector_particle) == radius * radius
 *       等价于函数体中定义的: a*(x^2) - 2b*x + c = 0
 * \param index           <[in] 颗粒索引
 * \param positionCenter  <[in] 颗粒中心
 * \param pointInside     <[in] 网格中心
 * \param pointOutside    <[in] 网格角点
 */
double mixVoidFraction::segmentParticleIntersection(int index,
                                                    vector positionCenter,
                                                    vector pointInside,
                                                    vector pointOutside) const {
  scalar lambda_1 = 0.0;
  scalar lambda_2 = 0.0;
  scalar lambda = 0.0;
  // 获取颗粒半径
  scalar radius = particleCloud_.radius(index);

  // 获取方程系数
  scalar a = (pointOutside - pointInside) & (pointOutside - pointInside);
  scalar c = (pointInside - positionCenter) & (pointInside - positionCenter);
  c -= radius * radius;
  scalar b = 2.0 * (pointOutside - pointInside) & (pointInside - positionCenter);
  double eps = 1.0e-12;

  scalar D = b * b - 4.0 * a * c;

  if (D >= 0.0) { // 如果方程有实数解, 则说明连线与球面一定有两个交点, 分别用 lambda_1 和 lambda_2 表示
    // 方程的两个实数解
    lambda_1 = (-1.0 * b + sqrt(D)) / (2.0 * a);
    lambda_2 = (-1.0 * b - sqrt(D)) / (2.0 * a);

    if (lambda_1 >= -eps && lambda_1 <= 1.0 + eps) {
      // 如果 lambda_1 在 [0, 1] 之间, 则说明这个交点在网格中心与网格角点的连线内, 则返回这个解
      lambda = lambda_1;
    }
    else if (lambda_2 >= -eps && lambda_2 <= 1.0 + eps) {
      // 如果 lambda_2 在 [0, 1] 之间, 则说明这个交点在网格中心与网格角点的连线内, 则返回这个解
      lambda = lambda_2;
    }
  }

  // 确保 lambda 在 [0, 1] 区间内, 而不是在 [-eps, 1 + eps] 区间内
  if (lambda > 1.0) {
    lambda = 1.0;
  }

  if (lambda < 0.0) {
    lambda = 0.0;
  }

  return lambda;
}

}  // End of namespace Foam
