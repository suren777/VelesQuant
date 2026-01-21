
void oneStepForwardADIDouglas(int t,
                              const std::vector<std::vector<double>> &inM,
                              std::vector<std::vector<double>> &outM) {
  oneStepForwardExplicit(t, inM, outM);
  std::vector<std::vector<double>> midM(outM);
  oneStepForwardDouglasX(t, inM, midM, outM);
  midM = outM;
  oneStepForwardDouglasY(t, inM, midM, outM);
}

void oneStepForwardExplicit(int t, const std::vector<std::vector<double>> &inM,
                            std::vector<std::vector<double>> &outM) {
  double te = gridT_[t + 1];
  double ti = gridT_[t];
  double delT = te - ti;
  double tm = 0.5 * (te + ti);
  double sigma1 = whichValue(tm, timeSigma1s_, sigma1s_);
  double difu1 = sigma1 * sigma1;
  double sigma2 = whichValue(tm, timeSigma2s_, sigma2s_);
  double difu2 = sigma2 * sigma2;
  double alpha = whichValue(tm, timeAlphas_, alphas_);

  // j=0
  double X0 = gridX_[0];
  double X1 = gridX_[1];
  double X2 = gridX_[2];
  double conv1 = -kappa1_ * X0;
  // j=0, k=0 lower boundary condition
  double Y0 = gridY_[0];
  double Y1 = gridY_[1];
  double Y2 = gridY_[2];
  double conv2 = -(lambda_ * X0 + kappa2_ * Y0);
  double rate = alpha + X0 + Y0;
  outM[0][0] =
      (1.0 / delT - rate) * inM[0][0] +
      (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
          inM[0][0] +
      (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][0] +
      (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][0] +
      (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
          inM[0][0] +
      (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[0][1] +
      (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[0][2];
  outM[0][0] = delT * outM[0][0];
  // j=0, k\=0 middle range
  double Y, Yu, Yl;
  for (int k = 1; k < Mxy_ - 1; k++) {
    Y = gridY_[k];
    Yu = gridY_[k + 1];
    Yl = gridY_[k - 1];
    conv2 = -(lambda_ * X0 + kappa2_ * Y);
    rate = alpha + X0 + Y;
    outM[0][k] =
        (1.0 / delT - rate) * inM[0][k] +
        (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
            inM[0][k] +
        (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][k] +
        (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][k] +
        (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[0][k - 1] +
        (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) * inM[0][k] +
        (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[0][k + 1];
    outM[0][k] = delT * outM[0][k];
  }
  // j=0, k=M-1 upper boundary condition
  double Y3l = gridY_[Mxy_ - 3];
  double Y2l = gridY_[Mxy_ - 2];
  double Y1l = gridY_[Mxy_ - 1];
  conv2 = -(lambda_ * X0 + kappa2_ * Y1l);
  rate = alpha + X0 + Y1l;
  outM[0][Mxy_ - 1] =
      (1.0 / delT - rate) * inM[0][Mxy_ - 1] +
      (-conv1 * (X2 + X1 - 2 * X0) + difu1) / (X1 - X0) / (X2 - X0) *
          inM[0][Mxy_ - 1] +
      (conv1 * (X2 - X0) - difu1) / (X1 - X0) / (X2 - X1) * inM[1][Mxy_ - 1] +
      (-conv1 * (X1 - X0) + difu1) / (X2 - X1) / (X2 - X0) * inM[2][Mxy_ - 1] +
      (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
          inM[0][Mxy_ - 3] +
      (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
          inM[0][Mxy_ - 2] +
      (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
          inM[0][Mxy_ - 1];
  outM[0][Mxy_ - 1] = delT * outM[0][Mxy_ - 1];

  // j\=0
  double X, Xu, Xl;
  for (int j = 1; j < Mxy_ - 1; j++) {
    X = gridX_[j];
    Xu = gridX_[j + 1];
    Xl = gridX_[j - 1];
    conv1 = -kappa1_ * X;
    // j\=0, k=0 lower boundary condition
    conv2 = -(lambda_ * X + kappa2_ * Y0);
    rate = alpha + X + Y0;
    outM[j][0] =
        (1.0 / delT - rate) * inM[j][0] +
        (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][0] +
        (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) * inM[j][0] +
        (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][0] +
        (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
            inM[j][0] +
        (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[j][1] +
        (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[j][2];
    outM[j][0] = delT * outM[j][0];
    // j\=0, k\=0 middle range
    for (int k = 1; k < Mxy_ - 1; k++) {
      Y = gridY_[k];
      Yu = gridY_[k + 1];
      Yl = gridY_[k - 1];
      conv2 = -(lambda_ * X + kappa2_ * Y);
      rate = alpha + X + Y;
      outM[j][k] =
          (1.0 / delT - rate) * inM[j][k] +
          (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) * inM[j - 1][k] +
          (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) / (Xu - X) *
              inM[j][k] +
          (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) * inM[j + 1][k] +
          (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) * inM[j][k - 1] +
          (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) / (Yu - Y) *
              inM[j][k] +
          (conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) * inM[j][k + 1];
      outM[j][k] = delT * outM[j][k];
    }
    // j\=0, k=M-1 upper boundary condition
    conv2 = -(lambda_ * X + kappa2_ * Y1l);
    rate = alpha + X + Y1l;
    outM[j][Mxy_ - 1] = (1.0 / delT - rate) * inM[j][Mxy_ - 1] +
                        (-conv1 * (Xu - X) + difu1) / (X - Xl) / (Xu - Xl) *
                            inM[j - 1][Mxy_ - 1] +
                        (conv1 * (Xu - 2 * X + Xl) - difu1) / (X - Xl) /
                            (Xu - X) * inM[j][Mxy_ - 1] +
                        (conv1 * (X - Xl) + difu1) / (Xu - X) / (Xu - Xl) *
                            inM[j + 1][Mxy_ - 1] +
                        (conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) /
                            (Y1l - Y3l) * inM[j][Mxy_ - 3] +
                        (-conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) /
                            (Y1l - Y2l) * inM[j][Mxy_ - 2] +
                        (conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) /
                            (Y1l - Y3l) * inM[j][Mxy_ - 1];
    outM[j][Mxy_ - 1] = delT * outM[j][Mxy_ - 1];
  }

  // j=M-1
  double X3l = gridX_[Mxy_ - 3];
  double X2l = gridX_[Mxy_ - 2];
  double X1l = gridX_[Mxy_ - 1];
  conv1 = -kappa1_ * X1l;
  // j=M-1, k=0 lower boundary condition
  conv2 = -(lambda_ * X1l + kappa2_ * Y0);
  rate = alpha + X1l + Y0;
  outM[Mxy_ - 1][0] =
      (1.0 / delT - rate) * inM[Mxy_ - 1][0] +
      (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
          inM[Mxy_ - 3][0] +
      (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
          inM[Mxy_ - 2][0] +
      (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
          inM[Mxy_ - 1][0] +
      (-conv2 * (Y2 + Y1 - 2 * Y0) + difu2) / (Y1 - Y0) / (Y2 - Y0) *
          inM[Mxy_ - 1][0] +
      (conv2 * (Y2 - Y0) - difu2) / (Y1 - Y0) / (Y2 - Y1) * inM[Mxy_ - 1][1] +
      (-conv2 * (Y1 - Y0) + difu2) / (Y2 - Y1) / (Y2 - Y0) * inM[Mxy_ - 1][2];
  outM[Mxy_ - 1][0] = delT * outM[Mxy_ - 1][0];
  // j=M-1, k\=0 middle range
  for (int k = 1; k < Mxy_ - 1; k++) {
    Y = gridY_[k];
    Yu = gridY_[k + 1];
    Yl = gridY_[k - 1];
    conv2 = -(lambda_ * X1l + kappa2_ * Y);
    rate = alpha + X1l + Y;
    outM[Mxy_ - 1][k] = (1.0 / delT - rate) * inM[Mxy_ - 1][k] +
                        (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) /
                            (X1l - X3l) * inM[Mxy_ - 3][k] +
                        (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) /
                            (X1l - X2l) * inM[Mxy_ - 2][k] +
                        (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) /
                            (X1l - X3l) * inM[Mxy_ - 1][k] +
                        (-conv2 * (Yu - Y) + difu2) / (Y - Yl) / (Yu - Yl) *
                            inM[Mxy_ - 1][k - 1] +
                        (conv2 * (Yu - 2 * Y + Yl) - difu2) / (Y - Yl) /
                            (Yu - Y) * inM[Mxy_ - 1][k] +
                        (-conv2 * (Y - Yl) + difu2) / (Yu - Y) / (Yu - Yl) *
                            inM[Mxy_ - 1][k + 1];
    outM[Mxy_ - 1][k] = delT * outM[Mxy_ - 1][k];
  }
  // j=M-1, k=M-1 upper boundary condition
  conv2 = -(lambda_ * X1l + kappa2_ * Y1l);
  rate = alpha + X1l + Y1l;
  outM[Mxy_ - 1][Mxy_ - 1] =
      (1.0 / delT - rate) * inM[Mxy_ - 1][Mxy_ - 1] +
      (-conv1 * (X1l - X2l) + difu1) / (X2l - X3l) / (X1l - X3l) *
          inM[Mxy_ - 3][Mxy_ - 1] +
      (conv1 * (X1l - X3l) - difu1) / (X2l - X3l) / (X1l - X2l) *
          inM[Mxy_ - 2][Mxy_ - 1] +
      (-conv1 * (2 * X1l - X2l - X3l) + difu1) / (X1l - X2l) / (X1l - X3l) *
          inM[Mxy_ - 1][Mxy_ - 1] +
      (-conv2 * (Y1l - Y2l) + difu2) / (Y2l - Y3l) / (Y1l - Y3l) *
          inM[Mxy_ - 1][Mxy_ - 3] +
      (conv2 * (Y1l - Y3l) - difu2) / (Y2l - Y3l) / (Y1l - Y2l) *
          inM[Mxy_ - 1][Mxy_ - 2] +
      (-conv2 * (2 * Y1l - Y2l - Y3l) + difu2) / (Y1l - Y2l) / (Y1l - Y3l) *
          inM[Mxy_ - 1][Mxy_ - 1];
  outM[Mxy_ - 1][Mxy_ - 1] = delT * outM[Mxy_ - 1][Mxy_ - 1];
}
