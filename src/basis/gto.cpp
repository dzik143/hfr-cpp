/******************************************************************************/
/*                                                                            */
/* Copyright (C) 2009 Sylwester Wysocki <sw143@wp.pl>                         */
/*                                                                            */
/* Source code available at: https://github.com/dzik143/hfr-cpp               */
/*                                                                            */
/* This program is free software: you can redistribute it and/or modify       */
/* it under the terms of the GNU General Public License as published by       */
/* the Free Software Foundation, either version 3 of the License, or          */
/* (at your option) any later version.                                        */
/*                                                                            */
/* This program is distributed in the hope that it will be useful,            */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of             */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              */
/* GNU General Public License for more details.                               */
/*                                                                            */
/* You should have received a copy of the GNU General Public License          */
/* along with this program. If not, see <http://www.gnu.org/licenses/>        */
/*                                                                            */
/******************************************************************************/

#include <error.h>
#include "gto.h"

using std::clog;
using std::endl;

//
// Load coefficients of GTO basis sets from the file.
//
// fname - file name to load (IN)
//

void BasisGTO::loadBasis(std::string &fname)
{
  // ===================
  // = Open input file =
  // ===================

  std::ifstream f(fname.c_str());

  if (!f.good())
  {
    throw string("error: can not to open file: ") + fname;
  }

  // ===============
  // = Log message =
  // ===============

  clog.fill('.');
  clog.width(35);
  clog << std::left << "Loading basis set";

  // ================================
  // = Load basis set from the file =
  // ================================

  while (!f.eof())
  {
    int Z;
    string sym;

    // Read Z-number (number of protons in nuclei).
    f >> sym;
    Z = Geometry::symToZ(sym);

    // =========================
    // = Next chemical element =
    // =========================

    if (sym[0] != '*')
    {
      // * is a comment
      int ilePowlok;
      f >> ilePowlok; // how many shells?

      while (ilePowlok--)
      {
        // =========
        // = Shell =
        // =========

        string typ; // s,p,d....

        int numberOfPrimitives; // how many function within contraction

        f >> typ >> numberOfPrimitives;

        if (typ == "SP")
        {
          // =====================================================
          // = Common SP contractions                             =
          // = The same parameters are used for S and P functions =
          // ======================================================

          // Create two gauss contractions.
          GTO s("S");
          GTO p("P");

          // Common coefficient both for S and P functions.
          double c;

          // Potentially different alpha for S and P functions.
          double as;
          double ap;

          // Load primitive GTOs within contraction.
          while (numberOfPrimitives--)
          {
            // GTO_S(c,a,r) = c *     exp(ar)
            // GTO_P(c,a,r) = c * r * exp(ar)
            f >> c >> as >> ap;
            s.addPrimitive(c, as);
            p.addPrimitive(c, ap);
          }

          // Apply functions to the element with given Z-number.
          kontr[Z].push_back(s);
          kontr[Z].push_back(p);
        }
        else
        {
          // ================
          // = General case =
          // ================

          // Create one contracted gauss.
          GTO psi(typ);

          double c, alfa;

          // Load primitives within contraction.
          while (numberOfPrimitives--)
          {
            f >> c >> alfa;
            psi.addPrimitive(c, alfa);
          }

          // Apply function to the element with given Z-number.
          kontr[Z].push_back(psi);
        }
      }
    }
  }

  clog << "O.K!" << endl;
  f.close();
}

//
// Dump basis set to human-readable text.
//

string BasisGTO::toString()
{
  std::ostringstream ss;

  vector<GTO>::iterator it;

  for (int Z = 0; Z < 100; Z++)
  {
    if (!kontr[Z].empty())
    {
      ss << Geometry::ZToSym(Z);

      for (it = kontr[Z].begin(); it != kontr[Z].end(); it++)
      {
        ss << "  " << *it << endl;
      }

      ss << endl;
    }
  }

  return ss.str();
}

//
// Create atomic orbitals set for concrete molecule.
//
// geom - molecule geometry (IN)
//

void BasisGTO::initAO(Geometry &geom)
{
  int ileAtom = geom.ileAtom();

  vector<GTO>::iterator it;

  clog.fill('.');
  clog.width(35);
  clog << std::left << "Generating atomic orbitals";

  // zliczaj prym. gaussy
  ilePrymGTO = 0;

  for (int i = 0; i < ileAtom; i++)
  {
    Atom atom = geom[i];
    int Z = atom.getZ();
    Vec3 atomR = atom.getR();

    if (Z > 2)
    {
      throw Error("only s integrals are implemented", "BASIS", "BasisGTO::initAO");
    }

    for (it = kontr[Z].begin(); it != kontr[Z].end(); it++)
    {
      // =============================
      // = Fetch data from basis set =
      // =============================

      int lTotal = it->getL();
      Vec3 psiR = it->getR();
      Vec3 R = atomR + psiR;
      vector<double> psiC = it->getC();
      vector<double> psiA = it->getA();

      ilePrymGTO += psiC.size();

      // ==========================
      // = Generate 2l+1 orbitals =
      // ==========================

      for (int lx = 0; lx <= lTotal; lx++)
      {
        for (int ly = 0; ly <= lTotal; ly++)
        {
          for (int lz = 0; lz <= lTotal; lz++)
          {
            if ((lx + ly + lz) <= lTotal)
            {
              AO.push_back(GTO(psiC, psiA, lx, ly, lz, R));
            }
          }
        }
      }
    }
  }

  int one = ilePrymGTO * ilePrymGTO;
  int two = int(pow(ilePrymGTO, 4.0) / 4.0);

  double mem = (one + two) * 4.0 / 1024 / 1024;

  // wypisuj 3 miejsa po przcinku.
  clog.precision(3);

  clog << "O.K!" << endl
       << endl
       << "----------------------------------------------------" << endl
       << "   Number of contracted GTO functions: " << AO.size() << endl
       << "   Number of primitive GTO functions:  " << ilePrymGTO << endl
       << "   Number of one-electron integrals:   " << (ilePrymGTO * ilePrymGTO) << endl
       << "   Number of two-electron integrals:   " << int(pow(ilePrymGTO, 4.0) / 4) << endl
       << "   Memory needed to store integrals:   " << mem << " MB" << endl
       << "----------------------------------------------------" << endl
       << endl;

}

BasisGTO::BasisGTO()
{
};

string BasisGTO::showAO()
{
}

int BasisGTO::getDimension()
{
  return AO.size();
}

//
// Calculate overlap matrix Spq = <p|q> for two contracted GTO functions.
//
// S - overlap matrix (OUT)
//

void BasisGTO::computeOverlap(Matrix &S)
{
  int i, j;

  double suma;
  double X;    // zm. pomocnicza
  double Rpq2; // odl. miedzy centrami

  vector<double> cp, cq;
  vector<double> wp, wq;

  int p, q;
  int lp, lq;
  int ileOrb;

  ileOrb = AO.size(); // ile orbitali

  for (p = 0; p < ileOrb; p++)
  {
    for (q = p; q < ileOrb; q++)
    {
      lp = AO[p].getL(); // pobierz typ
      lq = AO[q].getL(); // orbitalu (l)

      // =====================================
      // = Distance between function centers =
      // =====================================

      Vec3 rp, rq;         //
      rp = AO[p].getR();   // centra orbitali
      rq = AO[q].getR();   //
      Rpq2 = rp.dist2(rq); // kwadrat odleglosci

      // ==============================
      // = Calculate overlap integral =
      // ==============================

      // Sum all contractions
      suma = 0.0;

      cp = AO[p].getC(); // wspolczynniki
      cq = AO[q].getC(); // kontrakcji

      wp = AO[p].getA(); // wykladniki
      wq = AO[q].getA(); //

      for (i = 0; i < cp.size(); i++)
      {
        for (j = 0; j < cq.size(); j++)
        {
          X = (2.0 * wp[i] * wq[j]) / (wp[i] + wq[j]);

          suma = suma + cp[i] * cq[j]
                              * pow(wp[i] * wq[j], -0.75)
                              * pow(X, 1.5)
                              * exp(-0.5 * X * Rpq2);
        }
      }

      S(p, q) = suma;
      S(q, p) = suma;
    }
  }
}

//
// Calculate one-electron matrix:
// hpq = <p|h|q> = <p|T+V|q>
//
// h    - one-electron matrix (OUT)
// geom - molecule geometry (IN)
//

void BasisGTO::computeOneElectron(Matrix &h, Geometry &geom)
{
  int i, j;

  double suma;

  vector<double> cp, cq;
  vector<double> wp, wq;

  double X;             // zm. pomocnicza
  double Rpq, Rpq2, R2; // odl. miedzy centrami
  double dx, dy, dz;    // roznice wsp. cetrow
  double xa, ya, za;    // odl. od atomu
  double rox, roy, roz;
  double px, py, pz, qx, qy, qz;

  int ileOrb;
  int p, q;
  int lp, lq;
  int np, nq;
  int Ip, Iq;
  int Ilep, Ileq;
  int a; // numer atomu
  double Tpq, Vpq, S;
  double t; // pomocnicza do F0
  double c;

  ileOrb = AO.size();

  for (p = 0; p < ileOrb; p++)
  {
    for (q = p; q < ileOrb; q++)
    {

      lp = AO[p].getL(); // pobierz typ
      lq = AO[q].getL(); // orbitalu (l)

      // ===========================
      // = odleglosc miedzy centr. =
      // ===========================

      Vec3 rp, rq;         //
      rp = AO[p].getR();   // centra orbitali
      rq = AO[q].getR();   //
      Rpq2 = rp.dist2(rq); // kwadrat odleglosci

      // ==============
      // = licz calki =
      // ==============

      Tpq = 0.0;
      Vpq = 0.0;

      cp = AO[p].getC(); // wspolczynniki
      cq = AO[q].getC(); //   kontrakcji
      wp = AO[p].getA(); //
      wq = AO[q].getA(); // wykladniki

      for (i = 0; i < cp.size(); i++)
      {
        for (j = 0; j < cq.size(); j++)
        {

          // =========================
          // = licz calke nakrywania =
          // =========================

          X = (2.0 * wp[i] * wq[j]) / (wp[i] + wq[j]);

          S = cp[i] * cq[j] * pow(wp[i] * wq[j], -0.75)
                            * pow(X, 1.5)
                            * exp(-0.5 * X * Rpq2);

          // =========================
          // = licz calke en. kinet. =
          // =========================

          Tpq = Tpq + X * (1.5 - 0.5 * X * Rpq2) * S;

          // =========================
          // = licz calke prz. jadr. =
          // =========================

          px = AO[p].x();
          py = AO[p].y();
          pz = AO[p].z();

          qx = AO[q].x();
          qy = AO[q].y();
          qz = AO[q].z();

          //      wp*Rp + wq*Rq
          // ro = --------------
          //         wp + wq

          rox = (wp[i] * px + wq[j] * qx) / (wp[i] + wq[j]);
          roy = (wp[i] * py + wq[j] * qy) / (wp[i] + wq[j]);
          roz = (wp[i] * pz + wq[j] * qz) / (wp[i] + wq[j]);

          for (a = 0; a < geom.ileAtom(); a++)
          {
            xa = geom[a].x() - rox;
            ya = geom[a].y() - roy;
            za = geom[a].z() - roz;
            R2 = (xa * xa + ya * ya + za * za);

            t = (wp[i] + wq[j]) * R2;
            c = -1.12837916709551257 * geom[a].getZ() * sqrt(wp[i] + wq[j]);

            Vpq = Vpq + c * F0(t) * S;
          }
        }
      }

      h(p, q) = Tpq + Vpq;
      h(q, p) = Tpq + Vpq;
    }
  }
}

//
// Calculate two-electron matrix for concrete geometry:
// Grpqs = <pq|G|rs> = <pq| 1/r12 | rs>
//
// G - two-electron matrix (OUT)
//

void BasisGTO::computeTwoElectron(Tensor4 &G)
{
  int i, j, k, l;

  double suma;

  vector<double> wr_, ws_, wp_, wq_;
  vector<double> cr_, cs_, cp_, cq_;

  double wr, ws, wp, wq; // wykladniki p i q
  double cr, cs, cp, cq; // wsp. kontrakcji
  double X;              // zm. pomocnicza
  double Rpq2, Rrs2;     // odl. miedzy centrami
  double Ro2;            // odl. miedyz nowymi c.
  double dx, dy, dz;     // roznice wsp. cetrow
  double xrs, yrs, zrs;  // wsp. nowego centrum
  double xpq, ypq, zpq;  // wsp. nowego centrum
  double Rpq, Rrs;

  int lr, ls, lp, lq;
  int nr, ns, np, nq;
  int Ir, Is, Ip, Iq;
  int p, q, r, s;

  int ileOrb;

  double Srs, Spq;
  double t; // pomocnicza do F0
  double c;
  const double PI = 3.14159265358979324;

  ileOrb = AO.size();

  for (r = 0; r < ileOrb; r++)
  {
    for (s = r; s < ileOrb; s++)
    {
      for (p = 0; p < ileOrb; p++)
      {
        for (q = p; q < ileOrb; q++)
        {

          lr = AO[r].getL(); // pobierz typ
          ls = AO[s].getL(); // orbitalu (l)
          lp = AO[p].getL(); // orbitalu (l)
          lq = AO[q].getL(); // orbitalu (l)

          // ===========================
          // = licz odl. miedzy centr. =
          // ===========================

          Vec3 rr, rp, rq, rs;
          rr = AO[r].getR();
          rs = AO[s].getR();
          rp = AO[p].getR();
          rq = AO[q].getR();

          Rpq2 = rp.dist2(rq);
          Rrs2 = rr.dist2(rs);

          // ==============
          // = licz calki =
          // ==============

          wr_ = AO[r].getA();
          ws_ = AO[s].getA();
          wp_ = AO[p].getA();
          wq_ = AO[q].getA();

          cr_ = AO[r].getC();
          cs_ = AO[s].getC();
          cp_ = AO[p].getC();
          cq_ = AO[q].getC();

          suma = 0.0;

          for (i = 0; i < cr_.size(); i++)
          {
            wr = wr_[i];
            cr = cr_[i];

            for (j = 0; j < cs_.size(); j++)
            {
              ws = ws_[j];
              cs = cs_[j];

              for (k = 0; k < cp_.size(); k++)
              {
                wp = wp_[k];
                cp = cp_[k];

                for (l = 0; l < cq_.size(); l++)
                {
                  wq = wq_[l];
                  cq = cq_[l];

                  // =========================
                  // = licz calki nakrywania =
                  // =========================

                  X = (2.0 * wp * wq) / (wp + wq);

                  Spq = cp * cq * pow(wp * wq, -0.75)
                                * pow(X, 1.5) * exp(-0.5 * X * Rpq2);

                  X = (2.0 * wr * ws) / (wr + ws);

                  Srs = cr * cs * pow(wr * ws, -0.75)
                                * pow(X, 1.5) * exp(-0.5 * X * Rrs2);

                  // =========================
                  // = licz calke dwuelektr. =
                  // =========================

                  xrs = (wr * rr.x + ws * rs.x) / (wr + ws);
                  yrs = (wr * rr.y + ws * rs.y) / (wr + ws);
                  zrs = (wr * rr.z + ws * rs.z) / (wr + ws);

                  xpq = (wp * rp.x + wq * rq.x) / (wp + wq);
                  ypq = (wp * rp.y + wq * rq.y) / (wp + wq);
                  zpq = (wp * rp.z + wq * rq.z) / (wp + wq);

                  Ro2 = (xrs - xpq) * (xrs - xpq) +
                        (yrs - ypq) * (yrs - ypq) +
                        (zrs - zpq) * (zrs - zpq);

                  X = (wr + ws) * (wp + wq) / (wr + ws + wp + wq);
                  c = 1.12837916709551257 * sqrt(X); // 2*sqr(X/pi)
                  t = X * Ro2;

                  suma = suma + c * F0(t) * Srs * Spq;
                }
              }
            }
          }

          G(p, r, q, s) = suma;
          G(p, s, q, r) = suma;
          G(q, r, p, s) = suma;
          G(q, s, p, r) = suma;
        }
      }
    }
  }
}
