/* PROGRAM: tidecor
   VERSION: 2.5
   DATE :   December, 2010
   Author : Jason Chaffey


   Version 1.0 
 *******************************************************************************
 * This program for tide correction of hydrographic data is derived from
 * interp2d.c - Modifications and additions made by Patrick Roussel May 21,1999
 *****************************************************************************  

   Version 2.0   
   The FORTRAN components of Ver1.0 were translated to c for portability.
   Original FORTRAN written by M. Foreman at IOS.
   RAYBOUND bug fixed.
   Jason Chaffey              December 15, 2000
  

   Version 2.1   
   Added configuration file (tidecor2.?.cfg) so that the following are not
   hard-wired into the program:
       - mesh filenames (.nod and .ele);
       - the number of constituents to be used; and
       - for each constituent:
             - name of constituent; and
             - filename for data of that constituent.
   Jason Chaffey              January 15, 2001
  

   Version 2.2   
   Added support for complex input data (.v2c files)
   Jason Chaffey              March 23, 2001
  

   Version 2.3   
  
 * Tested against Pawlawicz's T_Tide and 
 * Foreman's tide4_r2  - Shawn Oakey
 *
  

   Version 2.4   
  
 * Debugging (float - double and int - long int conflicts) done
 * Switched to Manhattan distance vs "real" distance in closestnode
 *          - speed increase consideration
 * Thanks to Herman Varma for the above
 * Jason Chaffey              August 7, 2003
 *
  

   Version 2.42   
  
 * Bug in Manhattan distance calculation fixed. Returned to using it
 * for optimization.
 * Added fix for elements that cross International Dateline or
 * Greenwich Meridian in raybound.
 * Jason Chaffey              January 14, 2004
 *
  
   Version 2.5   
  
 * Fixed some problems with raybound for elements that cross International
 * Dateline or Greenwich Meridian in raybound. Problems became evident in the
 * global data set.
 * Jason Chaffey              June 16, 2008
 *
  
   Version 2.5.1   
  
 * Fixed a bug in raybound where triangle nodes were number 1-3 instead of
 * 0-2.
 * Jason Chaffey              November 29, 2010
 *
  
  
 * Dll Windows Version produced by Charles LeBlanc (Environment Canada Atlantic)
 * 
 * Charles LeBlanc              August, 2011
 *
*/

using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;

public partial class TidesAndCurrents
{
    #region Global Variables
    public DateTime StartDate { get; set; }
    public DateTime EndDate { get; set; }
    public double Longitude { get; set; }
    public double Latitude { get; set; }
    public double Steps { get; set; }
    public string DataPath { get; set; }
    public bool WLTrueCurFalse { get; set; }
    public bool Cancel { get; set; }

    public List<InputAndResult> InputAndResultList { get; set; }
    public List<Node> NodeList { get; set; }
    public List<Element> ElementList { get; set; }
    public List<Constituent> ConstituentNameList { get; set; }
    public List<MainConstituent> MainConstituentList { get; set; }
    public List<ShallowConstituent> ShallowConstituentList { get; set; }

    // Holder for Main Constituents
    public class MainConstituent
    {
        public MainConstituent()
        {
            dood = new List<int>();
            satellites = new List<Satellite>();
        }
        public string Name { get; set; }
        public List<int> dood { get; set; }
        public double phase { get; set; }
        public List<Satellite> satellites { get; set; }
        public double f { get; set; }
        public double freq { get; set; }
        public double vu { get; set; }

        // Holder for Satellites for each main constituent 
        public class Satellite
        {
            public Satellite()
            {
                deldList = new List<int>();
            }
            public List<int> deldList { get; set; }
            public double phase { get; set; }
            public double ratio { get; set; }
            public int corr { get; set; }
        }
    }
    // Holder for Shallow water constituents 
    public class ShallowConstituent
    {
        public ShallowConstituent()
        {
            factors = new List<Factor>();
        }
        public string Name { get; set; }
        public List<Factor> factors { get; set; }
        public double f { get; set; }
        public double freq { get; set; }
        public double vu { get; set; }

        // Holder for Constituent factors for the shallow water constituents 
        public class Factor
        {
            public Factor()
            {

            }
            public string Name { get; set; }
            public double factor { get; set; }
        }

    }
    // Holder for Nodes 
    public class Node
    {
        public Node()
        {

        }
        public int ID { get; set; }
        public double x { get; set; }
        public double y { get; set; }
        public double MinimumTide { get; set; }
        public double MinimumTide2 { get; set; }
    }
    // Holder for Elements
    public class Element
    {
        public Element()
        {
            this.NodeList = new List<Node>();
        }
        public int ID { get; set; }
        public List<Node> NodeList { get; set; }
    }
    // Holder for Constituent names and amp and phase values
    public class Constituent
    {
        public Constituent()
        {
            amp = new List<double>();
            phase = new List<double>();
            amp2 = new List<double>();
            phase2 = new List<double>();
        }
        public string Name { get; set; }
        public List<double> amp { get; set; }
        public List<double> phase { get; set; }
        public List<double> amp2 { get; set; }
        public List<double> phase2 { get; set; }
    }
    // Holder for Config file information
    private class Config
    {
        public Config()
        {

        }
        public string NodeFileName { get; set; }
        public string ElementFileName { get; set; }
        public string IOS_FileName { get; set; }
    }
    // Holder for location and date (memory input)
    public class InputAndResult
    {
        public InputAndResult()
        {

        }
        public double Longitude { get; set; }
        public double Latitude { get; set; }
        public DateTime Date { get; set; }
        public double Reslt { get; set; }
        public double Reslt2 { get; set; }
    }
    // Holder for Data Files name and Path
    private class DataFilePath
    {
        public DataFilePath()
        {

        }
        public string Name { get; set; }
        public string Path { get; set; }
    }
    // Holder for Astro information and parameters
    private class Astro
    {
        public Astro()
        {

        }
        public double d1 { get; set; }
        public double h { get; set; }
        public double pp { get; set; }
        public double s { get; set; }
        public double p { get; set; }
        public double enp { get; set; }
        public double dh { get; set; }
        public double dpp { get; set; }
        public double ds { get; set; }
        public double dp { get; set; }
        public double dnp { get; set; }
    }

    List<double> basis;                   // Factors for interpolation 
    Node closestNode = new Node();        // Closest node number to an arbitrary point 
    Config config = new Config();
    List<string> ParsedValues = new List<string>();

    #endregion Global Variables

    #region Constructors
    public TidesAndCurrents(DateTime startDate, DateTime endDate, double longitude, double latitude, 
        double steps, string dataPath, bool wLTrueCurFalse)
    {
        StartDate = startDate;
        EndDate = endDate;
        Longitude = longitude;
        Latitude = latitude;
        Steps = steps;
        DataPath = dataPath;
        WLTrueCurFalse = wLTrueCurFalse;
    }
    #endregion Constructors

    #region Functions
    private MainConstituent.Satellite AddSatellite(string satin)
    {
        /* Gets satellite info from input line, puts the info into a new
            satellite structure and returns the structure. */
        MainConstituent.Satellite newsat = new MainConstituent.Satellite();

        ParsedValues = satin.Split(default(string[]), StringSplitOptions.RemoveEmptyEntries).ToList<string>();

        newsat.deldList.Add(int.Parse(ParsedValues[0]));
        newsat.deldList.Add(int.Parse(ParsedValues[1]));
        newsat.deldList.Add(int.Parse(ParsedValues[2]));
        newsat.phase = double.Parse(ParsedValues[3]);
        if (ParsedValues[4].Contains("R"))
        {
            newsat.ratio = double.Parse(ParsedValues[4].Substring(0, ParsedValues[4].IndexOf("R")));
        }
        else
        {
            newsat.ratio = double.Parse(ParsedValues[4]);
        }

        if (!satin.Contains("R"))
        {
            newsat.corr = 0;
        }
        else
        {
            if (satin.Contains("R1"))
            {
                newsat.corr = 1;
            }
            else if (satin.Contains("R2"))
            {
                newsat.corr = 2;
            }
            else
            {
                newsat.corr = 0;
            }
        }
        return newsat;
    }
    private void AstroAngles(ref Astro astro)
    {
        /* Calculates the following ephermides of the sun and moon:
            astro.h  = mean longitude of the sun;
            astro.pp = mean longitude of the solar perigee;
            astro.s  = mean longitude of the moon;
            astro.p  = mean longitude of the lunar perigee; and
            astro.enp = negative of the longitude of the mean ascending node.
            Also calculates their rate of change.
            Units are cycles ( cycles / 365 days for rates of change ).
            These formulae were taken from pp 98 and 107 of the Explanatory
            Supplement to the Astronomical Ephermeris and the American
            Ephermis and Nautical Almanac (1961) */
        double d12, d2, d22, d23, f, f2;

        d12 = astro.d1 * astro.d1;
        d2 = astro.d1 * 1.0E-04;
        d22 = d2 * d2;
        d23 = Math.Pow(d2, 3.0);
        f = 360.0;
        f2 = f / 365.0;

        astro.h = (2.79696678E+02 + astro.d1 * 9.856473354E-01 + d22 * 2.267E-05) / f;
        astro.h -= (int)astro.h;

        astro.pp = (2.81220833E+02 + astro.d1 * 4.70684E-05 + d22 * 3.39E-05 + d23 * 7.0E-08) / f;
        astro.pp -= (int)astro.pp;

        astro.s = (2.70434164E+02 + astro.d1 * 1.31763965268E+01 - d22 * 8.5E-05 + d23 * 3.9E-08) / f;
        astro.s -= (int)astro.s;

        astro.p = (3.34329556E+02 + astro.d1 * 1.114040803E-01 - d22 * 7.739E-04 - d23 * 2.6E-07) / f;
        astro.p -= (int)astro.p;

        astro.enp = (-2.59183275E+02 + astro.d1 * 5.29539222E-02 - d22 * 1.557E-04 - d23 * 5.0E-08) / f;
        astro.enp -= (int)astro.enp;

        astro.dh = (9.856473354E-01 + astro.d1 * 2.267E-05 * 2.0E-08) / f2;

        astro.dpp = (4.70684E-05 + astro.d1 * 3.39E-05 * 2.0E-08 + d12 * 7.0E-08 * 3.0E-12) / f2;

        astro.ds = (1.31763965268E+01 - astro.d1 * 8.5E-05 * 2.0E-08 + d12 * 3.9E-08 * 3.0E-12) / f2;

        astro.dp = (1.114040803E-01 - astro.d1 * 7.739E-04 * 2.0E-08 - d12 * 2.6E-07 * 3.0E-12) / f2;

        astro.dnp = (5.29539222E-02 - astro.d1 * 1.557E-04 * 2.0E-08 - d12 * 5.0E-08 * 3.0E-12) / f2;
    }
    private Element Basis2d(InputAndResult mi)
    {
        basis = new List<double>();

        // Finds the closest node to a point (ptx, pty) and the element containing
        // that point, if one exists.
        // Also gets the basis functions for interpolations to that point.
        // Returns the element number contaning the point (ptx, pty).
        // Returns -999 if no element found containing the point (ptx, pty ). 

        int flag;
        List<double> xlocal = new List<double>();
        List<double> ylocal = new List<double>();

        for (int i = 0; i < 3; i++)
        {
            xlocal.Add(0.0);
            ylocal.Add(0.0);
        }

        closestNode = FindClosestNode(mi);

        //  Try to find an element that contains the point and the closest node.
        Element ele = new Element();
        ele = null;

        foreach (Element element in ElementList.Where(el => el.NodeList[0] == closestNode
                                        || el.NodeList[1] == closestNode
                                        || el.NodeList[2] == closestNode).Distinct<Element>().ToList<Element>())
        {
            int i = 0;
            foreach (Node n in element.NodeList)
            {
                xlocal[i] = n.x;
                ylocal[i] = n.y;
                i += 1;
            }

            // See if the point is within this element 
            flag = RayBound(xlocal, ylocal, mi);
            if (flag == 1)
            {
                // The point is within the element 
                ele = element;
                basis = Phi2d(xlocal, ylocal, mi);
                break;
            }
        }

        /* if the closest node's elements don't work, search through all elements */
        if (ele == null)
        {
            foreach (Element element in ElementList)
            {
                int i = 0;
                foreach (Node n in element.NodeList)
                {
                    xlocal[i] = n.x;
                    ylocal[i] = n.y;
                    i += 1;
                }

                /* See if the point is within this element */
                flag = RayBound(xlocal, ylocal, mi);
                if (flag == 1)
                {
                    ele = element;
                    basis = Phi2d(xlocal, ylocal, mi);
                    return (element);
                }
            }
        }
        return (ele);
    }
    private bool CalculateResults()
    {
        DateTime CurrentDate = DateTime.Now;
        Element ele = new Element();
        ele = null;
        List<double> elem_res = new List<double>();
        List<double> elem_res2 = new List<double>();
        List<double> elem_min = new List<double>();
        List<double> elem_min2 = new List<double>();

        // Loop through the input file. For each line calculate the tidal correction 

        RaiseMessageEvent(string.Format("Calculating Results for {0} date and time", InputAndResultList.Count));
        RaiseStatusEvent(string.Format("Calculating Results for {0} date and time", InputAndResultList.Count));

        foreach (InputAndResult mi in InputAndResultList)
        {
            if (Cancel)
            {
                Cancel = false;
                break;
            }

            if (CurrentDate.Day != mi.Date.Day)
            {
                RaiseStatusEvent(string.Format("{0:yyyy/MM/dd}", mi.Date));
            }

            ele = Basis2d(mi);

            if (ele == null)
            {
                // No element was found that contained the new position. 
                // Calculate the tidal correction for the closest node. 

                // removed by Charles LeBlanc, results should be -999 if no elements were found
                mi.Reslt = -999;
                mi.Reslt2 = -999;

                RaiseErrorEvent(string.Format("Wrong data path [{0}]. Please verify you data source.", DataPath));
                return false;
                //__________________________________________________
                //reslt = TideP(mi, closestNode, true);
                //reslt = reslt + closestNode.MinimumTide;
                //if (dataType == DataType.Current)
                //{
                //    reslt2 = TideP(mi, closestNode, false);
                //    reslt2 = reslt2 + closestNode.MinimumTide2;
                //}
                //__________________________________________________

            }
            else
            {
                /* If the point is inside and element, calculate the tidal correction for
                    each node of the element and interpolate to get the tidal correction
                    at the new position. */
                elem_res.Clear();
                elem_res2.Clear();
                elem_min.Clear();
                elem_min2.Clear();

                foreach (Node n in ele.NodeList)
                {
                    elem_res.Add(TideP(mi, n, true));
                    elem_min.Add(n.MinimumTide);
                    if (!WLTrueCurFalse)
                    {
                        elem_res2.Add(TideP(mi, n, false));
                        elem_min2.Add(n.MinimumTide2);
                    }
                }

                mi.Reslt = elem_res[0] * basis[0]
                        + elem_res[1] * basis[1]
                    + elem_res[2] * basis[2]
                    + elem_min[0] * basis[0]
                        + elem_min[1] * basis[1]
                    + elem_min[2] * basis[2];
                if (!WLTrueCurFalse)
                {
                    mi.Reslt2 = elem_res2[0] * basis[0]
                            + elem_res2[1] * basis[1]
                        + elem_res2[2] * basis[2]
                        + elem_min2[0] * basis[0]
                            + elem_min2[1] * basis[1]
                        + elem_min2[2] * basis[2];
                }
            }
        }

        return true;
    }
    private bool CheckAllDataOK()
    {
        if (Longitude < -180 || Longitude > 180)
        {
            RaiseErrorEvent("Please enter a valid Longitude. Needs to be between -180 and 180.");
            return false;
        }

        if (Latitude < -90 || Latitude > 90)
        {
            RaiseErrorEvent("Please enter a valid Latitude. Needs to be between -90 and 90.");
            return false;
        }

        if (Steps < 1)
        {
            RaiseErrorEvent("Please enter a valid Steps value. Needs to be >= 1");
            return false;
        }

        if (StartDate.Year < 1900)
        {
            RaiseErrorEvent("Please enter a valid Start Date. Needs to be > 1900");
            return false;
        }

        if (StartDate >= EndDate)
        {
            RaiseErrorEvent("Please enter an end date > than the start date.");
            return false;
        }

        return true;
    }
    private bool CreateMemoryInput()
    {
        InputAndResultList = new List<InputAndResult>();

        RaiseMessageEvent("Creating Memory Input ...");
        RaiseStatusEvent("Creating Memory Input ...");

        while (EndDate > StartDate)
        {
            InputAndResultList.Add(new InputAndResult()
            {
                Longitude = Longitude,
                Latitude = Latitude,
                Date = StartDate
            });
            StartDate = StartDate.AddMinutes(Steps);
        }

        RaiseMessageEvent(string.Format("{0} Memory Input Created", InputAndResultList.Count));

        return true;
    }
    private Node FindClosestNode(InputAndResult mi)
    {
        // Find the node that is closest to the point (mi.ptx, mi.pty). 
        double currdist, closedist;
        Node n = new Node();
        n = null;

        closedist = Math.Abs(mi.Latitude - NodeList[0].y) + Math.Abs(mi.Longitude - NodeList[0].x);

        foreach (Node node in NodeList)
        {
            currdist = Math.Abs(mi.Latitude - node.y) + Math.Abs(mi.Longitude - node.x);
            if (currdist < closedist)
            {
                closedist = currdist;
                n = node;
            }
        }
        return (n);
    }
    private long GetJulianDay(DateTime Date)
    {
        /* Calculate the Julian day number.
            Accounts for the change to the Gregorian calandar. */
        long jul, ja, jy, jm;

        jy = Date.Year;
        if (jy == 0)
        {
            RaiseErrorEvent("JULDAY: There is no year 0!");
            return 0;
        }
        if (jy < 0) ++jy;
        if (Date.Month > 2)
        {
            jm = Date.Month + 1;
        }
        else
        {
            --jy;
            jm = Date.Month + 13;
        }

        jul = (long)(Math.Floor(365.25 * jy) + Math.Floor(30.6001 * jm) +
                        Date.Day + 1720995);
        if ((Date.Day + 31L * (Date.Month + 12L * Date.Year)) >= (15 + (long)31 * (10 + (long)12 * 1582)))
        {
            ja = (long)(0.01 * jy);
            jul += 2 - ja + (long)(0.25 * ja);
        }
        return jul;
    }
    public bool LoadBase()
    {
        if (!ReadConfig())
            return false;

        if (!ReadNodes())
            return false;

        if (!ReadElements())
            return false;

        if (!ReadConstituents())
            return false;

        if (!ReadIOS_TideTable())
            return false;

        if (!LoadConstituentsValues())
            return false;

        return true;
    }
    private bool LoadConstituentsValues()
    {
        foreach (Constituent ct in ConstituentNameList)
        {
            StreamReader sr;

            FileInfo fiConstituent = new FileInfo(DataPath + ct.Name + ".barotropic." + (WLTrueCurFalse ? "s2c" : "v2c"));

            RaiseMessageEvent(string.Format("Reading Constituent Values [{0}]", fiConstituent.FullName));
            RaiseStatusEvent(string.Format("Reading Constituent Values [{0}]", fiConstituent.FullName));
            
            if (!fiConstituent.Exists)
            {
                RaiseErrorEvent(string.Format("Constituent file does not exist [{0}].", fiConstituent.FullName));
                return false;
            }

            RaiseStatusEvent(string.Format("Reading file [{0}].", fiConstituent.FullName));

            try
            {
                sr = fiConstituent.OpenText();
            }
            catch (Exception ex)
            {
                RaiseMessageEvent(string.Format("Error message [{0}].", ex.Message));
                RaiseErrorEvent(string.Format("Error while trying to open [{0}].", fiConstituent.FullName));
                return false;
            }

            sr.ReadLine();
            sr.ReadLine();
            if (ct.Name != "Z0")
            {
                sr.ReadLine();
            }

            foreach (Node node in NodeList)
            {
                ParsedValues = sr.ReadLine().Split(default(string[]), StringSplitOptions.RemoveEmptyEntries).ToList<string>();
                if (WLTrueCurFalse)
                {
                    if (ParsedValues.Count() != 3)
                    {
                        RaiseErrorEvent(string.Format("[{0}] file not parsed properly.", fiConstituent.FullName));
                        return false;
                    }
                    ct.amp.Add(double.Parse(ParsedValues[1]));
                    ct.phase.Add(double.Parse(ParsedValues[2]));
                }
                else
                {
                    if (ct.Name == "Z0")
                    {
                        if (ParsedValues.Count() != 3)
                        {
                            RaiseErrorEvent(string.Format("[{0}] file not parsed properly.", fiConstituent.FullName));
                            return false;
                        }
                        ct.amp.Add(double.Parse(ParsedValues[1]));
                        ct.amp2.Add(double.Parse(ParsedValues[2]));
                        ct.phase.Add((double)0.0);
                        ct.phase2.Add((double)0.0);
                    }
                    else
                    {
                        if (ParsedValues.Count() != 5)
                        {
                            RaiseErrorEvent(string.Format("[{0}] file not parsed properly.", fiConstituent.FullName));
                            return false;
                        }
                        ct.amp.Add(double.Parse(ParsedValues[1]));
                        ct.phase.Add(double.Parse(ParsedValues[2]));
                        ct.amp2.Add(double.Parse(ParsedValues[3]));
                        ct.phase2.Add(double.Parse(ParsedValues[4]));
                    }
                }
                if (sr.EndOfStream)
                {
                    break;
                }
            }

            sr.Close();
        }
        return true;
    }
    private List<double> Phi2d(List<double> xloc, List<double> yloc, InputAndResult mi)
    {
        /* Calculates the basis functions for interpolating to a point inside
            an element. */
        List<double> phi = new List<double>();
        double area, a, b, c;
        int j = 0, k = 0;

        area = 0.5 * (xloc[0] * (yloc[1] - yloc[2]) +
                        xloc[1] * (yloc[2] - yloc[0]) +
                        xloc[2] * (yloc[0] - yloc[1]));

        /* Calculate the Basis function... */
        for (int i = 0; i <= 2; i++)
        {
            switch (i)
            {
                case 0: j = 1; k = 2; break;
                case 1: j = 2; k = 0; break;
                case 2: j = 0; k = 1; break;
            }
            a = (xloc[j] * yloc[k] - xloc[k] * yloc[j]) / (area * 2);
            b = (yloc[j] - yloc[k]) / (area * 2);
            c = -1 * (xloc[j] - xloc[k]) / (area * 2);
            phi.Add(a + b * mi.Longitude + c * mi.Latitude);
        }
        return (phi);
    }
    private int RayBound(List<double> xd, List<double> yd, InputAndResult mi)
    {
        /*  Subroutine to check wether or not a point is inside a polygon.
        The process is as follows:
            Use an arbitrary ray (here, y = constant and x >= xref), starting from 
        the point and going off to infinity.
            Count the number of polygon boundaries it crosses.
            If an odd number, the point is inside the polygon, otherwise it is
        outside.   */
        int j, bcross;
        double b, m, x;

        bcross = 0; /* Number of boundary crossings. */

        /* Check to see if the element side crosses the International Dateline
            (changes sign at +180/-180 degrees) and if so, change the longitudes
            so that they all have the same sign. */

        if (mi.Longitude > 0.0)
        {
            if ((xd[0] < -170.0) && ((xd[1] > 170.0) || (xd[2] > 170.0)))
            {
                xd[0] += 360.0;
            }
            if ((xd[1] < -170.0) && ((xd[0] > 170.0) || (xd[2] > 170.0)))
            {
                xd[1] += 360.0;
            }
            if ((xd[2] < -170.0) && ((xd[1] > 170.0) || (xd[0] > 170.0)))
            {
                xd[2] += 360.0;
            }
        }
        else
        {
            if ((xd[0] > 170.0) && ((xd[1] < -170.0) || (xd[2] < -170.0)))
            {
                xd[0] -= 360.0;
            }
            if ((xd[1] > 170.0) && ((xd[0] < -170.0) || (xd[2] < -170.0)))
            {
                xd[1] -= 360.0;
            }
            if ((xd[2] > 170.0) && ((xd[1] < -170.0) || (xd[0] < -170.0)))
            {
                xd[2] -= 360.0;
            }
        }

        /* As above, except for the Greenwich meridian, for longitude coordinates
            that range from 0 to 360 degrees. */

        if (mi.Longitude > 350.0)
        {
            if ((xd[0] < 10.0) && ((xd[1] > 350.0) || (xd[2] > 350.0)))
            {
                xd[0] += 360.0;
            }
            if ((xd[1] < 10.0) && ((xd[0] > 350.0) || (xd[2] > 350.0)))
            {
                xd[1] += 360.0;
            }
            if ((xd[2] < 10.0) && ((xd[1] > 350.0) || (xd[0] > 350.0)))
            {
                xd[2] += 360.0;
            }
        }
        else
        {
            if ((xd[0] > 350.0) && ((xd[1] < 10.0) || (xd[2] < 10.0)))
            {
                xd[0] -= 360.0;
            }
            if ((xd[1] > 350.0) && ((xd[0] < 10.0) || (xd[2] < 10.0)))
            {
                xd[1] -= 360.0;
            }
            if ((xd[2] > 350.0) && ((xd[1] < 10.0) || (xd[0] < 10.0)))
            {
                xd[2] -= 360.0;
            }
        }

        for (int i = 0; i <= 2; i++)
        {

            /* for each line segment around the element */
            j = ((i == 2) ? 0 : i + 1);

            /* If both endpoints of the line segment are on the same (vertical)
            side of the ray, do nothing.
                Otherwise, count the number of times the ray intersects the segment. */
            if (!(((yd[i] < mi.Latitude) && (yd[j] < mi.Latitude)) ||
                    ((yd[i] >= mi.Latitude) && (yd[j] >= mi.Latitude))))
            {

                if (xd[i] != xd[j])
                {
                    m = (yd[j] - yd[i]) / (xd[j] - xd[i]);
                    b = yd[i] - m * xd[i];
                    x = (mi.Latitude - b) / m;
                    if (x > mi.Longitude)
                    {
                        bcross++;
                    }
                }
                else
                {
                    if (xd[j] > mi.Longitude)
                    {
                        bcross++;
                    }
                }
            }
        }

        /*  Return the evenness/oddness of the boundary crossings
                i.e. the remainder from division by two. */
        return (bcross % 2);
    }
    private bool ReadConfig()
    {
        string TheLine = "";
        string Last4Letters = "";
        StreamReader sr;

        FileInfo fiCfg = new FileInfo(DataPath + "tidecor.cfg");

        RaiseMessageEvent(string.Format("Reading Config File [{0}]", fiCfg.FullName));
        RaiseStatusEvent(string.Format("Reading Config File [{0}]", fiCfg.FullName));

        if (!fiCfg.Exists)
        {
            RaiseErrorEvent(string.Format("Config file could not be found [{0}].", fiCfg.FullName));
            return false;
        }

        try
        {
            sr = fiCfg.OpenText();
        }
        catch (Exception ex)
        {
            RaiseMessageEvent(string.Format("Error message {0}.", ex.Message));
            RaiseErrorEvent(string.Format("Could not open file {0}.", fiCfg.FullName));
            return false;
        }

        while (true)
        {
            TheLine = sr.ReadLine().Trim();
            Last4Letters = TheLine.Substring(TheLine.Length - 4).ToLower();
            switch (Last4Letters)
            {
                case ".nod":
                    {
                        config.NodeFileName = TheLine;
                    }
                    break;
                case ".ele":
                    {
                        config.ElementFileName = TheLine;
                    }
                    break;
                case "etbl":
                    {
                        config.IOS_FileName = TheLine;
                    }
                    break;
                default:
                    break;
            }
            if (sr.EndOfStream)
            {
                break;
            }
        }

        sr.Close();

        return true;
    }
    private bool ReadConstituents()
    {
        ConstituentNameList = new List<Constituent>();
        StreamReader sr;

        FileInfo fiConst = new FileInfo(DataPath + "constituents.txt");

        RaiseMessageEvent(string.Format("Reading Constituent Names [{0}]",fiConst.FullName));
        RaiseStatusEvent(string.Format("Reading Constituent Names [{0}]", fiConst.FullName));

        if (!fiConst.Exists)
        {
            RaiseErrorEvent(string.Format("Could not find the constituents file [{0}].", fiConst.FullName));
            return false;
        }

        try
        {
            sr = fiConst.OpenText();
        }
        catch (Exception ex)
        {
            RaiseMessageEvent(string.Format("Error message [{0}].", ex.Message));
            RaiseErrorEvent(string.Format("Error while trying to open [{0}].", fiConst.FullName));
            return false;
        }

        while (true)
        {
            ParsedValues = sr.ReadLine().Split(default(string[]), StringSplitOptions.RemoveEmptyEntries).ToList<string>();
            if (ParsedValues.Count() != 1)
            {
                RaiseErrorEvent("Constituents file not parsed properly.");
                return false;
            }
            if (ParsedValues[0].Trim().ToUpper() == "NONE")
            {
                break;
            }
            ConstituentNameList.Add(new Constituent() { Name = ParsedValues[0].Trim().ToUpper() });
            if (sr.EndOfStream)
            {
                break;
            }
        }
        sr.Close();
        return true;
    }
    private bool ReadElements()
    {
        ElementList = new List<Element>();
        StreamReader sr;

        FileInfo fiElem = new FileInfo(DataPath + config.ElementFileName);

        RaiseMessageEvent(string.Format("Reading Element File [{0}]", fiElem.FullName));
        RaiseStatusEvent(string.Format("Reading Element File [{0}]", fiElem.FullName));

        if (!fiElem.Exists)
        {
            RaiseErrorEvent(string.Format("Could not find the element file [{0}].", fiElem.FullName));
            return false;
        }

        try
        {
            sr = fiElem.OpenText();
        }
        catch (Exception ex)
        {
            RaiseMessageEvent(string.Format("Error message [{0}].", ex.Message));
            RaiseErrorEvent(string.Format("Error while trying to open [{0}].", fiElem.FullName));
            return false;
        }

        while (true)
        {
            ParsedValues = sr.ReadLine().Split(default(string[]), StringSplitOptions.RemoveEmptyEntries).ToList<string>();
            if (ParsedValues.Count() != 4)
            {
                RaiseErrorEvent("Element file not parsed properly.");
                return false;
            }
            List<Node> TempNodeList = new List<Node>();
            TempNodeList.Add(NodeList[int.Parse(ParsedValues[1]) - 1]);
            TempNodeList.Add(NodeList[int.Parse(ParsedValues[2]) - 1]);
            TempNodeList.Add(NodeList[int.Parse(ParsedValues[3]) - 1]);

            ElementList.Add(new Element() { ID = int.Parse(ParsedValues[0]), NodeList = TempNodeList });
            if (sr.EndOfStream)
            {
                break;
            }
        }
        sr.Close();

        return true;
    }
    private bool ReadIOS_TideTable()
    {
        MainConstituentList = new List<MainConstituent>();
        ShallowConstituentList = new List<ShallowConstituent>();

        /* Read in the constituent data from the IOS_tidetbl file */
        int cnt, nln;
        string satin1, satin2, satin3;

        FileInfo fiIOS = new FileInfo(DataPath + config.IOS_FileName);

        RaiseMessageEvent(string.Format("Reading IOS Tide Table [{0}]", fiIOS.FullName));
        RaiseStatusEvent(string.Format("Reading IOS Tide Table [{0}]", fiIOS.FullName));

        if (!fiIOS.Exists)
        {
            RaiseErrorEvent(string.Format("Could not find file [{0}].", fiIOS.FullName));
            return false;
        }

        /* Read in the main constituents*/

        StreamReader sr = fiIOS.OpenText();

        while (true) // this while will not hit the end of stream, it will be stopped by a blank line.
        {
            ParsedValues = sr.ReadLine().Split(default(string[]), StringSplitOptions.RemoveEmptyEntries).ToList<string>();

            if (ParsedValues.Count == 0)
            {
                break;
            }
            MainConstituent TempNewCon = new MainConstituent();
            TempNewCon.Name = ParsedValues[0];
            List<int> dood = new List<int>();
            TempNewCon.dood.Add(int.Parse(ParsedValues[1]));
            TempNewCon.dood.Add(int.Parse(ParsedValues[2]));
            TempNewCon.dood.Add(int.Parse(ParsedValues[3]));
            TempNewCon.dood.Add(int.Parse(ParsedValues[4]));
            TempNewCon.dood.Add(int.Parse(ParsedValues[5]));
            if (ParsedValues[6].LastIndexOf("-") > 0)
            {
                ParsedValues.Add(ParsedValues[7]);
                ParsedValues[7] = ParsedValues[6].Substring(ParsedValues[6].LastIndexOf("-"));
                ParsedValues[6] = ParsedValues[6].Substring(0, ParsedValues[6].LastIndexOf("-"));
            }
            TempNewCon.dood.Add(int.Parse(ParsedValues[6]));
            TempNewCon.phase = double.Parse(ParsedValues[7]);

            int NumberOfSatellite = int.Parse(ParsedValues[8]);
            /* Read in the satellites for this constituent, if any */
            if (NumberOfSatellite > 0)
            {
                nln = ((NumberOfSatellite - 1) / 3) + 1;
                for (int i = 0; i < nln; i++)
                {
                    string TheLine = sr.ReadLine();
                    cnt = NumberOfSatellite - (i * 3); /* # of satellites on this line */
                    switch (cnt)
                    {
                        case 1:
                            satin1 = TheLine.Substring(12);
                            TempNewCon.satellites.Add(AddSatellite(satin1));
                            break;
                        case 2:
                            satin1 = TheLine.Substring(12, 23);
                            satin2 = TheLine.Substring(12 + 23);
                            TempNewCon.satellites.Add(AddSatellite(satin1));
                            TempNewCon.satellites.Add(AddSatellite(satin2));
                            break;
                        default:
                            satin1 = TheLine.Substring(12, 23);
                            satin2 = TheLine.Substring(12 + 23, 23);
                            satin3 = TheLine.Substring(12 + 23 + 23);
                            TempNewCon.satellites.Add(AddSatellite(satin1));
                            TempNewCon.satellites.Add(AddSatellite(satin2));
                            TempNewCon.satellites.Add(AddSatellite(satin3));
                            break;
                    }
                }
            }
            else
            {
                /*  NO satellites */
                TempNewCon.satellites = null;
            }
            MainConstituentList.Add(TempNewCon);

            if (sr.EndOfStream)
            {
                break;
            }
            if (ParsedValues.Count() == 0)
            {
                return true;
            }
        }



        /* Read in the shallow water constiuents */

        while (true)
        {
            ParsedValues = sr.ReadLine().Split(default(string[]), StringSplitOptions.RemoveEmptyEntries).ToList<string>();

            if (ParsedValues.Count == 0)
            {
                break;
            }

            ShallowConstituent ns = new ShallowConstituent();

            ns.Name = ParsedValues[0];
            int NumberOfFactors = int.Parse(ParsedValues[1]);

            switch (NumberOfFactors)
            {
                case 4:
                    {
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[3], factor = double.Parse(ParsedValues[2]) });
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[5], factor = double.Parse(ParsedValues[4]) });
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[7], factor = double.Parse(ParsedValues[6]) });
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[9], factor = double.Parse(ParsedValues[8]) });
                    }
                    break;
                case 3:
                    {
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[3], factor = double.Parse(ParsedValues[2]) });
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[5], factor = double.Parse(ParsedValues[4]) });
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[7], factor = double.Parse(ParsedValues[6]) });
                    }
                    break;
                case 2:
                    {
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[3], factor = double.Parse(ParsedValues[2]) });
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[5], factor = double.Parse(ParsedValues[4]) });
                    }
                    break;
                case 1:
                    {
                        ns.factors.Add(new ShallowConstituent.Factor() { Name = ParsedValues[3], factor = double.Parse(ParsedValues[2]) });
                    }
                    break;
            }

            ShallowConstituentList.Add(ns);

            if (sr.EndOfStream)
            {
                break;
            }
        }

        sr.Close();

        return true;
    }
    private bool ReadNodes()
    {
        NodeList = new List<Node>();
        StreamReader sr;
        FileInfo fiNode = new FileInfo(DataPath + config.NodeFileName);

        RaiseMessageEvent(string.Format("Reading Node File [{0}]", fiNode.FullName));
        RaiseStatusEvent(string.Format("Reading Node File [{0}]", fiNode.FullName));

        if (!fiNode.Exists)
        {
            RaiseErrorEvent(string.Format("Could not find the node file [{0}].", fiNode.FullName));
            return false;
        }

        try
        {
            sr = fiNode.OpenText();
        }
        catch (Exception ex)
        {
            RaiseMessageEvent(string.Format("Error message [{0}].", ex.Message));
            RaiseErrorEvent(string.Format("Error while trying to open [{0}].", fiNode.FullName));
            return false;
        }

        while (true)
        {
            ParsedValues = sr.ReadLine().Split(default(string[]), StringSplitOptions.RemoveEmptyEntries).ToList<string>();
            if (ParsedValues.Count() != 3)
            {
                RaiseErrorEvent("Node file not parsed properly.");
                return false;
            }
            NodeList.Add(new Node() { ID = int.Parse(ParsedValues[0]), x = double.Parse(ParsedValues[1]), y = double.Parse(ParsedValues[2]), MinimumTide = 0, MinimumTide2 = 0 });
            if (sr.EndOfStream)
            {
                break;
            }
        }

        sr.Close();

        return true;
    }
    public bool RunMain(bool NeedToReload)
    {
        if (!CheckAllDataOK())
            return false;

        if (NeedToReload)
        {
            if (!LoadBase())
                return false;
        }

        if (!CreateMemoryInput())
            return false;

        if (!CalculateResults())
            return false;

        RaiseMessageEvent("Successfuly completed ...");
        RaiseStatusEvent("Successfuly completed ...");
        return true;
    }
    private void SetVuf(long kh, double xlat)
    {
        Astro astro = new Astro();

        /* Calculate the amplitudes, phases, etc. for each of the constituents */
        double hh, tau, dtau, slat, sumc, sums, v, vdbl;
        double adj = 0, uu, uudbl;
        long kd, ktmp;

        kd = GetJulianDay(new DateTime(1899, 12, 31));
        astro.d1 = (double)(kh - kd) - 0.5;
        ktmp = kh * 24;


        AstroAngles(ref astro);

        hh = (double)ktmp - (Math.Floor((double)(ktmp / 24.0)) * 24.0);
        tau = hh / 24.0 + astro.h - astro.s;
        dtau = 365.00 + astro.dh - astro.ds;
        slat = Math.Sin((Math.PI) * (xlat / 180.0));

        /* The main constituents */
        foreach (MainConstituent mc in MainConstituentList)
        {
            mc.freq = (mc.dood[0] * dtau + mc.dood[1] * astro.ds + mc.dood[2] * astro.dh + mc.dood[3] * astro.dp +
                                mc.dood[4] * astro.dnp + mc.dood[5] * astro.dpp) / (24.0 * 365.0);

            vdbl = mc.dood[0] * tau + mc.dood[1] * astro.s + mc.dood[2] * astro.h + mc.dood[3] * astro.p +
                                mc.dood[4] * astro.enp + mc.dood[5] * astro.pp + mc.phase;

            v = vdbl - (Math.Floor(Math.Floor(vdbl) / 2.0) * 2.0);

            sumc = 1.0;
            sums = 0.0;

            if (mc.satellites != null)
            {
                foreach (MainConstituent.Satellite satellite in mc.satellites)
                {
                    switch (satellite.corr)
                    {
                        case 0:
                            adj = satellite.ratio;
                            break;
                        case 1:
                            adj = satellite.ratio * 0.36309 * (1.0 - 5.0 * slat * slat) / slat;
                            break;
                        case 2:
                            adj = satellite.ratio * 2.59808 * slat;
                            break;
                    }
                    uudbl = (double)satellite.deldList[0] * astro.p + (double)satellite.deldList[1] * astro.enp +
                            (double)satellite.deldList[2] * astro.pp + (double)satellite.phase;
                    uu = uudbl - (int)uudbl;

                    sumc += (adj * Math.Cos(uu * Math.PI * 2));
                    sums += (adj * Math.Sin(uu * Math.PI * 2));
                }
            }
            mc.f = Math.Sqrt((sumc * sumc) + (sums * sums));
            mc.vu = v + Math.Atan2(sums, sumc) / (Math.PI * 2);
        }

        /* The shallow water constituents */
        foreach (ShallowConstituent shallowConstituent in ShallowConstituentList)
        {
            shallowConstituent.f = 1.0;
            shallowConstituent.vu = 0.0;
            shallowConstituent.freq = 0.0;
            foreach (ShallowConstituent.Factor factor in shallowConstituent.factors)
            {
                foreach (MainConstituent mainConstituent in MainConstituentList)
                {
                    if (mainConstituent.Name == factor.Name)
                    {
                        shallowConstituent.f *= Math.Pow(mainConstituent.f, Math.Abs(factor.factor));
                        shallowConstituent.vu += (factor.factor * mainConstituent.vu);
                        shallowConstituent.freq += (factor.factor * mainConstituent.freq);
                        break;
                    }
                }
            }
        }
    }
    private double TideP(InputAndResult mi, Node n, bool IsWL)
    {
        /* Calculates and returns the tidal correction */
        int indx;
        long kd;
        double dthr, radgmt, revgmt, res;

        kd = GetJulianDay(mi.Date);
        SetVuf(kd, mi.Latitude);
        dthr = (((double)mi.Date.Hour * 3600.0) + ((double)mi.Date.Minute * 60.0) +
                    (double)mi.Date.Second) / 3600.0;
        res = 0.0;

        /* For each of the desired constituents ... (See top of program) */
        foreach (Constituent constituent in ConstituentNameList)
        {
            /* Find the constituent from those loaded from IOS_tidetbl */
            indx = Vuf(constituent.Name);
            if (indx < 0)
            {
                RaiseErrorEvent(string.Format("Bad Input Constituent: {0}.", constituent));
                return -999;
            }

            if (indx < MainConstituentList.Count)
            {                        // Main constituent 
                if (IsWL)
                    revgmt = MainConstituentList[indx].freq * dthr + MainConstituentList[indx].vu - constituent.phase[n.ID - 1] / 360.0;
                else
                    revgmt = MainConstituentList[indx].freq * dthr + MainConstituentList[indx].vu - constituent.phase2[n.ID - 1] / 360.0;

                radgmt = Math.PI * 2 * (revgmt - (int)revgmt);

                if (IsWL)
                    res += MainConstituentList[indx].f * constituent.amp[n.ID - 1] * Math.Cos(radgmt);
                else
                    res += MainConstituentList[indx].f * constituent.amp2[n.ID - 1] * Math.Cos(radgmt);
            }
            else if ((indx - MainConstituentList.Count) < ShallowConstituentList.Count)
            {     // Shallow water constituent 
                indx -= MainConstituentList.Count;
                if (IsWL)
                    revgmt = ShallowConstituentList[indx].freq * dthr + ShallowConstituentList[indx].vu - constituent.phase[n.ID - 1] / 360.0;
                else
                    revgmt = ShallowConstituentList[indx].freq * dthr + ShallowConstituentList[indx].vu - constituent.phase2[n.ID - 1] / 360.0;

                radgmt = Math.PI * 2 * (revgmt - (int)revgmt);
                if (IsWL)
                    res += ShallowConstituentList[indx].f * constituent.amp[n.ID - 1] * Math.Cos(radgmt);
                else
                    res += ShallowConstituentList[indx].f * constituent.amp2[n.ID - 1] * Math.Cos(radgmt);
            }
            else
            {
                RaiseErrorEvent(string.Format("Error in index: {0} {1} {2}.", indx, MainConstituentList.Count, ShallowConstituentList.Count));
                return -999;
            }
        }

        return (res);
    }
    private int Vuf(string inname)
    {
        /* Finds constituent info corresponding to inname and returns the index to
            the node containing the info. Shallow water constituent indices are
            returned as their number greater than the max # of main constituents.
            e.g. if we want the 2nd shallow water constituent and there are 45
            main constituents, then the indice returned is 46, since constituents
            are counted from zero. ( 45 - 1 + 2 = 46 ) */
        int i, j;

        i = 0;
        j = 0;
        while (MainConstituentList[i].Name != inname)
        {
            i++;
            if (i == MainConstituentList.Count) break;
        }
        if (i == MainConstituentList.Count)
        {
            j = 0;
            while (ShallowConstituentList[j].Name != inname)
            {
                j++;
                if (j == ShallowConstituentList.Count) break;
            }
        }

        if (i < MainConstituentList.Count)
        {
            return (i);
        }
        else if (j < ShallowConstituentList.Count)
        {
            return (i + j);
        }
        else
        {
            RaiseErrorEvent(string.Format("Constituent %s not found!", inname));
            return -999;
        }
    }
    #endregion Functions

    // Miscellaneous messages event ex: beginning file read
    #region Message Event
    // event attributes
    public class MessageEventArgs : EventArgs
    {
        // Properties
        public string Message { get; set; }

        // Constructor
        public MessageEventArgs(string Message)
        {
            this.Message = Message;

        }
    }
    public class ErrorEventArgs : EventArgs
    {
        // Properties
        public string Error { get; set; }

        public ErrorEventArgs(string Error)
        {
            this.Error = Error;
        }
    }
    public class StatusEventArgs : EventArgs
    {
        // Properties
        public string Status { get; set; }

        public StatusEventArgs(string Status)
        {
            this.Status = Status;
        }
    }
    
    // delegate function attached to in parent object
    public delegate void TidesAndCurrentsMessageEventHandler(object sender, MessageEventArgs e);
    public delegate void TidesAndCurrentsErrorEventHandler(object sender, ErrorEventArgs e);
    public delegate void TidesAndCurrentsStatusEventHandler(object sender, StatusEventArgs e);

    // event function called from within TidesAndCurrents
    public event TidesAndCurrentsMessageEventHandler TidesAndCurrentsMessageEvent;
    public event TidesAndCurrentsErrorEventHandler TidesAndCurrentErrorEvent;
    public event TidesAndCurrentsStatusEventHandler TidesAndCurrentStatusEvent;

    private void RaiseMessageEvent(string Message)
    {
        if (this.TidesAndCurrentsMessageEvent != null)
            TidesAndCurrentsMessageEvent(this, new MessageEventArgs(Message));
    }
    private void RaiseErrorEvent(string Error)
    {
        if (this.TidesAndCurrentErrorEvent != null)
            TidesAndCurrentErrorEvent(this, new ErrorEventArgs(Error));
    }
    private void RaiseStatusEvent(string Status)
    {
        if (this.TidesAndCurrentStatusEvent != null)
            TidesAndCurrentStatusEvent(this, new StatusEventArgs(Status));
    }

    #endregion

}
