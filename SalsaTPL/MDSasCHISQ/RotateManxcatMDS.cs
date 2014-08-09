using System;
using System.IO;
using System.Text;
using System.Threading.Tasks;
using Manxcat;
using SALSALibrary;
using Salsa.Core.Configuration.Sections;

namespace RotateMDS
{
    public class RotateManxcatMDS
    {
        // Parameters in 3D case
        // 0 Translate-X
        // 1 Translate-Y
        // 2 Translate-Z
        // 3 THETA-X
        // 4 THETA-Y
        // 5 THETA-Z
        // 6 Scale (Negative if Inversion)
        //
        // Parameters in 2D case
        // 0 Translate-X
        // 1 Translate-Y
        // 2 THETA-Z
        // 3 Scale (Negative if Reflection in x-axis)

        public static SALSAFileProperties RotationFileProperties;
        public static SALSADataPointProperties[] RotationPointProperties;

        public static double[][] FirstData; // Initial point data
        public static double[][] SecondData; // Second point data
        public static double[] FirstMean; // Mean of initial point data
        public static double[] SecondMean; // Mean of initial point data

        public static int ScaleOption = 0;
        public static int RotationOption = 0;
        public static int RotationOption_invert = 0;
        public static int RotationOption_GenerateTestInput = 0;
        public static bool InvertSolution = false; // If true Invert solution
        public static double FirstScale; // Scale of first data
        public static double SecondScale; // Scale of first data
        public static int PointVectorDimension = 3; // Vector Dimension
        public static int ScalePosition = -1; // Position in Parameter Array of Scale
        public static double[,] CurrentCosmicScaling; // Current Matrix Multiplication for Inversion
        public static double[][,] SavedCosmicScaling; // Saved Matrix Multiplication for Inversion
        public static int NumberSubRotations = 3; // Number of 2D Rotation matrices
        public static double MinimumDistance = .001; // Minimum Distance for Scaling

        public static double ScaleA1 = 1.0;
        public static double ScaleA2 = 0.5;
        public static double ScaleAlpha = 0.75;


        public static void SetupHotsunforRotateMDS()
        {
            PointVectorDimension = ManxcatCentral.Configuration.LocalVectorDimension;
            if (PointVectorDimension == 3)
            {
                Hotsun.npar = 7;
            }
            else if (PointVectorDimension == 2)
            {
                Hotsun.npar = 4;
            }
            else
            {
                Exception e = SALSAUtility.SALSAError("Unsupported Dimension " + PointVectorDimension.ToString());
            }
            Hotsun.ndata = SALSAUtility.PointCount_Global*PointVectorDimension;
            Hotsun.ParameterVectorDimension = 1;
            GenericManxcat.SetupGenericManxcat();
            return;
        }

        // End SetupHotsunforMDS

        public static void SetupRotateMDS()
        {
            if (Hotsun.InitializationLoops < 2)
            {
                Hotsun.InitializationLoops = 2;
                Hotsun.InitLoopChisq = new double[Hotsun.InitializationLoops];
                ManxcatCentral.Configuration.InitializationLoops = Hotsun.InitializationLoops;
            }
            Hotsun.bnderrLimit = 10;

            if (SALSAUtility.NumberFixedPoints > 0)
            {
                Exception e = SALSAUtility.SALSAError(" Fixed Points not Allowed in Rotations");
                throw (e);
            }
            RotationOption = ManxcatCentral.Configuration.RotationOption; // set Rotation Option
            RotationOption_GenerateTestInput = RotationOption/100;
                // Control special test inputs = 0 normal -- don't generate points to rotate;
            //  RotationOption_GenerateTestInput = 1 full inversion; 
            //  RotationOption_GenerateTestInput = 2  invert x and y only; 
            //  RotationOption_GenerateTestInput = 3 invert x only
            //  RotationOption_GenerateTestInput = 4 set new file = old file
            RotationOption = RotationOption - 100*RotationOption_GenerateTestInput;
            RotationOption_invert = RotationOption/10;
                // Control how to control scaling in fit. = 0 force scaling positive so when cosmic scaling improper, one cannot fit proper rotation
            // RotationOption_invert = 1 allow scaling to be negative which can allow fit to be trapped in wring proper/improper domain for all starts. 
            // If RotationOption_invert=0 then every other initialization loop is proper (0 1 ..) or improper (1 3 ...)
            //  Note initial angle values are 0 for initialization loops 0 1nd 1, they are random for later loops
            //  Angles are forced to lie between -PI  and PI
            RotationOption = RotationOption - 10*RotationOption_invert;
            // RotationOption = 0 usual Chisq; = 1 Each term divided by magnitude of original point so one minimizes realative errors squared

            PointVectorDimension = SALSAUtility.GlobalFileProperties.LocalVectorDimension;
            if (PointVectorDimension != SALSAUtility.GlobalFileProperties.LocalVectorDimension)
            {
                Exception e =
                    SALSAUtility.SALSAError(" Inconsistent Small Dimension in Rotations " +
                                            PointVectorDimension.ToString());
                throw (e);
            }

            int InitializationFileType = -1;
            int InitializationNumberofPoints = -1;

            //  Set Reference Points
            if (SALSA_ProcessVariedandFixed.AreValuesSet(SALSAUtility.GlobalPointProperties))
            {
                SALSAUtility.SALSAPrint(2,
                                        " Reference Points Taken from basic Labels File " +
                                        ManxcatCentral.Configuration.DataLabelsFileName);
            }
            else
            {
                // Set up Values of Used Points from Initialization File
                var InitializationFileProperties = new SALSAFileProperties();
                var InitializationPointProperties = new SALSADataPointProperties[SALSAUtility.NumberOriginalPoints];

                string OriginalMDSFileName = ManxcatCentral.Configuration.InitializationFileName;
//                if (!OriginalMDSFileName.Contains(":") || !OriginalMDSFileName.Contains("$"))
//                    OriginalMDSFileName = ManxcatCentral.Configuration.ControlDirectoryName + "\\" + OriginalMDSFileName;
                SALSAUtility.SALSAPrint(2, " Reference Points Taken from Initialization File " + OriginalMDSFileName);

                // Begin Changes saliya 03/25/11
                // Note. I think we can agree on the zero based indices. The following should be removed then.
                if (OriginalMDSFileName.Contains("SIMPLE"))
                    InitializationFileProperties.LocalPointStartIndex = 1;
                // End Changes saliya 03/25/11

                SALSA_Properties.ReadDataPointFile(OriginalMDSFileName, ref InitializationFileType,
                                                   InitializationFileProperties, ref InitializationPointProperties,
                                                   ref InitializationNumberofPoints);
                if ((SALSAUtility.NumberOriginalPoints < InitializationNumberofPoints) ||
                    (InitializationFileProperties.NumberOriginalPoints != SALSAUtility.NumberOriginalPoints))
                {
                    Exception e =
                        SALSAUtility.SALSAError(" Inconsistent Initialization File Point Counts " +
                                                InitializationNumberofPoints.ToString() + " or "
                                                + InitializationFileProperties.NumberOriginalPoints.ToString() +
                                                " Expected is " + SALSAUtility.NumberOriginalPoints.ToString());
                    throw (e);
                }

                for (int InitialIndex = 0; InitialIndex < InitializationNumberofPoints; InitialIndex++)
                {
                    int OriginalIndex = InitializationPointProperties[InitialIndex].OriginalPointNumber;
                    if (!InitializationPointProperties[InitialIndex].valuesset)
                        continue;
                    SALSAUtility.GlobalPointProperties[OriginalIndex].x = InitializationPointProperties[InitialIndex].x;
                    SALSAUtility.GlobalPointProperties[OriginalIndex].y = InitializationPointProperties[InitialIndex].y;
                    SALSAUtility.GlobalPointProperties[OriginalIndex].z = InitializationPointProperties[InitialIndex].z;

                    SALSAUtility.GlobalPointProperties[OriginalIndex].valuesset =
                        InitializationPointProperties[InitialIndex].valuesset;
                }


                if (!SALSA_ProcessVariedandFixed.AreValuesSet(SALSAUtility.GlobalPointProperties))
                {
                    Exception e = SALSAUtility.SALSAError(" Reference Points Not set");
                    throw (e);
                }
            }

            // Now Read Points to be Rotated
            // Set up Values of Used Points from Rotation File
            int RotationFileType = -1;
            RotationFileProperties = new SALSAFileProperties();
            RotationPointProperties = new SALSADataPointProperties[SALSAUtility.NumberOriginalPoints];
            if (RotationOption_GenerateTestInput == 0)
            {
                int RotationNumberofPoints = -1;
                string RotationFileName = ManxcatCentral.Configuration.RotationLabelsFileName;
//                if (!RotationFileName.Contains(":"))
//                    RotationFileName = ManxcatCentral.Configuration.ControlDirectoryName + "\\" + RotationFileName;
                if (RotationFileName.Contains("SIMPLE"))
                    RotationFileProperties.LocalPointStartIndex = 1;
                SALSAUtility.SALSAPrint(2, " Points to be rotated Taken from File " + RotationFileName);
                SALSA_Properties.ReadDataPointFile(RotationFileName, ref RotationFileType, RotationFileProperties,
                                                   ref RotationPointProperties, ref RotationNumberofPoints);
                if ((SALSAUtility.NumberOriginalPoints < RotationNumberofPoints) ||
                    (RotationFileProperties.NumberOriginalPoints != SALSAUtility.NumberOriginalPoints))
                {
                    Exception e =
                        SALSAUtility.SALSAError(" Inconsistent Rotation File Point Counts " +
                                                RotationNumberofPoints.ToString() + " or "
                                                + RotationFileProperties.NumberOriginalPoints.ToString() +
                                                " Expected is " + SALSAUtility.NumberOriginalPoints.ToString());
                    throw (e);
                }
                if (!SALSA_ProcessVariedandFixed.AreValuesSet(RotationPointProperties))
                {
                    Exception e = SALSAUtility.SALSAError(" Points to rotate Not set");
                    throw (e);
                }
            } // End normal case when we read file of points to be rotated
            else
            {
                // Generate file to be Rotated
                RotationFileType = InitializationFileType;
                for (int PointIndex = 0; PointIndex < SALSAUtility.NumberOriginalPoints; PointIndex++)
                {
                    RotationPointProperties[PointIndex] = new SALSADataPointProperties();
                }

                double fudgex = 1.0;
                double fudgey = 1.0;
                double fudgez = 1.0;
                if (RotationOption_GenerateTestInput == 1)
                {
                    SALSAUtility.SALSAPrint(2, " Points to be rotated generated as -x -y -z ");
                    fudgex = -1.0;
                    fudgey = -1.0;
                    fudgez = -1.0;
                }
                if (RotationOption_GenerateTestInput == 2)
                {
                    SALSAUtility.SALSAPrint(2, " Points to be rotated generated as -x -y +z ");
                    fudgex = -1.0;
                    fudgey = -1.0;
                }
                if (RotationOption_GenerateTestInput == 3)
                {
                    SALSAUtility.SALSAPrint(2, " Points to be rotated generated as -x +y +z ");
                    fudgex = -1.0;
                }
                if (RotationOption_GenerateTestInput >= 4)
                {
                    SALSAUtility.SALSAPrint(2, " Points to be rotated generated as +x +y +z ");
                }
                for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
                {
                    int OriginalIndex = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
                    RotationPointProperties[OriginalIndex].x = fudgex*
                                                               SALSAUtility.GlobalPointProperties[OriginalIndex].x;
                    ;
                    RotationPointProperties[OriginalIndex].y = fudgey*
                                                               SALSAUtility.GlobalPointProperties[OriginalIndex].y;
                    RotationPointProperties[OriginalIndex].z = fudgez*
                                                               SALSAUtility.GlobalPointProperties[OriginalIndex].z;
                    RotationPointProperties[OriginalIndex].valuesset =
                        SALSAUtility.GlobalPointProperties[OriginalIndex].valuesset;
                }
            }

            //  Set up operational data
            FirstData = new double[SALSAUtility.PointCount_Global][]; // Initial point data
            SecondData = new double[SALSAUtility.PointCount_Global][];
            ; // Second point data
            FirstMean = new double[PointVectorDimension]; // Mean of initial point data
            SecondMean = new double[PointVectorDimension];
            for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
            {
                FirstData[GlobalPointIndex] = new double[PointVectorDimension];
                SecondData[GlobalPointIndex] = new double[PointVectorDimension];
            }

            //  Process Data
            ProcessPointData(SALSAUtility.GlobalPointProperties, FirstData, FirstMean, out FirstScale);
            ProcessPointData(RotationPointProperties, SecondData, SecondMean, out SecondScale);
            ScalePosition = 6;
            if (PointVectorDimension == 2)
                ScalePosition = 3;

            NumberSubRotations = 3;
            if (PointVectorDimension == 2)
                NumberSubRotations = 1;

            // This holds the inversion transformation
            CurrentCosmicScaling = new double[PointVectorDimension,PointVectorDimension];
            SavedCosmicScaling = new double[Hotsun.InitializationLoops][,];
            for (int scalingindex = 0; scalingindex < Hotsun.InitializationLoops; scalingindex++)
            {
                SavedCosmicScaling[scalingindex] = new double[PointVectorDimension,PointVectorDimension];
            }

            // If ManxcatCentral.MetadataforRun.MinimumDistance  is positive, it is absolute Minimum Distance
            // If ManxcatCentral.MetadataforRun.MinimumDistance  is negative, it is - multiplier of System Average
            if (ManxcatCentral.Configuration.MinimumDistance < 0)
                MinimumDistance = -ManxcatCentral.Configuration.MinimumDistance*FirstScale;
            else
                MinimumDistance = ManxcatCentral.Configuration.MinimumDistance;

            string firstM = "";
            string secondM = "";
            for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
            {
                firstM += FirstMean[LocalVectorIndex].ToString("E4") + " ";
                secondM += SecondMean[LocalVectorIndex].ToString("E4") + " ";
            }
            SALSAUtility.SALSAPrint(0,
                                    " First Scale " + FirstScale.ToString("E4") + " Second Scale " +
                                    SecondScale.ToString("E4") + " First Mean " + firstM + " Second Mean " + secondM);
        }

        // End SetupRotateMDS()

        public static bool Calcfg(Desertwind Solution)
        {
            // call Accum with value and derivative for each point

            bool violat = false;
            double Scale = 1.0;
            double DerivScale = 1.0;

            if (ScaleOption == 0)
            {
                if ((Solution.param[ScalePosition][0] < 0.0) && (RotationOption_invert == 0))
                {
                    violat = true;
                    Solution.param[ScalePosition][0] = -Solution.param[ScalePosition][0];
                }
                Solution.param[ScalePosition][0] = Math.Max(ScaleAlpha*FirstScale/SecondScale,
                                                            Solution.param[ScalePosition][0]);
                Scale = Solution.param[ScalePosition][0];
            }
            else
            {
                Solution.param[ScalePosition][0] = ProperAngle(Solution.param[ScalePosition][0]);
                Scale = CalculateScale(Solution.param[ScalePosition][0]);
                DerivScale = DifferentiateScale(Solution.param[ScalePosition][0]);
                // SALSAUtility.SALSAPrint(0, "Scale DerivScale " + Scale.ToString("E4") + " " + DerivScale.ToString("E4"));
            }

            //  Force Angles between -PI and PI
            if (PointVectorDimension == 3)
            {
                for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
                {
                    Solution.param[LocalVectorIndex + PointVectorDimension][0] =
                        ProperAngle(Solution.param[LocalVectorIndex + PointVectorDimension][0]);
                }
            }
            if (PointVectorDimension == 2)
            {
                Solution.param[2][0] = ProperAngle(Solution.param[2][0]);
            } // End Force Angles between -PI and PI

            //  Each Value is weight * (Scaling * Cosmic * Rotation * Second Vector + Translation - First Vector)
            //  Each Derivative is weight * ( Scaling * Cosmic * Rotation Term * Second Vectot + Treanslation Term )
            //  Rotation, Rotation Term, Translation, Translation Term are independent of Point
            //  weight depends on Option Used 

            var Translation = new double[PointVectorDimension];
            var Rotation = new double[PointVectorDimension,PointVectorDimension];
            var DerivTranslation = new double[Hotsun.npar][];
            var DerivRotationScale = new double[Hotsun.npar][,];
            var SubRotation = new double[NumberSubRotations][,];
            var SubDeriveRotation = new double[NumberSubRotations][,];

            for (int ipar = 0; ipar < Hotsun.npar; ipar++)
            {
                DerivTranslation[ipar] = new double[PointVectorDimension];
                DerivRotationScale[ipar] = new double[PointVectorDimension,PointVectorDimension];
                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
                {
                    DerivTranslation[ipar][LocalVectorIndex1] = 0.0;
                    for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                    {
                        if (ipar < PointVectorDimension)
                            DerivRotationScale[ipar][LocalVectorIndex1, LocalVectorIndex2] = 0.0;
                                // Rotation Matrix Zero for Translations
                    }
                }
            }

            for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
            {
                Translation[LocalVectorIndex] = Solution.param[LocalVectorIndex][0]; // Translations
                DerivTranslation[LocalVectorIndex][LocalVectorIndex] = 1.0; // Deriv of Translation
            }

            for (int SubRotationIndex = 0; SubRotationIndex < NumberSubRotations; SubRotationIndex++)
            {
                SubRotation[SubRotationIndex] = new double[PointVectorDimension,PointVectorDimension];
                SubDeriveRotation[SubRotationIndex] = new double[PointVectorDimension,PointVectorDimension];
                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
                {
                    for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                    {
                        SubDeriveRotation[SubRotationIndex][LocalVectorIndex1, LocalVectorIndex2] = 0.0;
                        SubRotation[SubRotationIndex][LocalVectorIndex1, LocalVectorIndex2] = 0.0;
                    }
                }
                SetRotationMatrix(SubRotation[SubRotationIndex], SubRotationIndex,
                                  Solution.param[PointVectorDimension + SubRotationIndex][0], false);
                SetRotationMatrix(SubDeriveRotation[SubRotationIndex], SubRotationIndex,
                                  Solution.param[PointVectorDimension + SubRotationIndex][0], true);
            }
            if (PointVectorDimension == 3)
            {
                SetTotalRotation(Rotation, Scale, CurrentCosmicScaling, SubRotation[2], SubRotation[1], SubRotation[0]);
                SetTotalRotation(DerivRotationScale[PointVectorDimension], Scale, CurrentCosmicScaling, SubRotation[2],
                                 SubRotation[1], SubDeriveRotation[0]); // Rotation for First Angle
                SetTotalRotation(DerivRotationScale[PointVectorDimension + 1], Scale, CurrentCosmicScaling,
                                 SubRotation[2], SubDeriveRotation[1], SubRotation[0]); // Rotation for Second Angle
                SetTotalRotation(DerivRotationScale[PointVectorDimension + 2], Scale, CurrentCosmicScaling,
                                 SubDeriveRotation[2], SubRotation[1], SubRotation[0]); // Rotation for Third Angle
                SetTotalRotation(DerivRotationScale[ScalePosition], DerivScale, CurrentCosmicScaling, SubRotation[2],
                                 SubRotation[1], SubRotation[0]); // Rotation for Scale
            }
            if (PointVectorDimension == 2)
            {
                SetTotalRotation(Rotation, Scale, CurrentCosmicScaling, SubRotation[0]);
                SetTotalRotation(DerivRotationScale[PointVectorDimension], Scale, CurrentCosmicScaling,
                                 SubDeriveRotation[0]); // Rotation for First Angle
                SetTotalRotation(DerivRotationScale[ScalePosition], DerivScale, CurrentCosmicScaling, SubRotation[0]);
                    // Rotation for Scale
            }

            //  Set up Accum
            GenericManxcat.ReInitializeAccum();

            //  Finally Loop over entries in Chisq
            Parallel.For(0, SALSAUtility.ParallelOptions.MaxDegreeOfParallelism, SALSAUtility.ParallelOptions,
                         (ThreadNo) =>
                             {
                                 var DerivativeFl = new double[PointVectorDimension][];
                                 for (int LocalVectorIndex = 0;
                                      LocalVectorIndex < PointVectorDimension;
                                      LocalVectorIndex++)
                                 {
                                     DerivativeFl[LocalVectorIndex] = new double[Hotsun.npar];
                                 }

                                 int indexlen = SALSAUtility.PointsperThread[ThreadNo];
                                 int beginpoint = SALSAUtility.StartPointperThread[ThreadNo] -
                                                  SALSAUtility.PointStart_Process;

                                 for (int DistributedPointIndex = beginpoint;
                                      DistributedPointIndex < indexlen + beginpoint;
                                      DistributedPointIndex++)
                                 {
                                     int GlobalIndex = SALSAUtility.PointStart_Process + DistributedPointIndex;

                                     // Setting weight to Math.Sqrt(1.0/SALSAUtility.PointCount_Global)* Fudge Factor/FirstScale
                                     double weight = Math.Sqrt(1.0/SALSAUtility.PointCount_Global)*
                                                     ManxcatCentral.ChisqFunctionCalcMultiplier/FirstScale;

                                     // RotationOption =0 is simple least squares sum; RotationOption = 1 implies use fractional errors
                                     if (RotationOption == 1)
                                     {
                                         double size = 0.0;
                                         for (int LocalVectorIndex = 0;
                                              LocalVectorIndex < PointVectorDimension;
                                              LocalVectorIndex++)
                                         {
                                             size += FirstData[GlobalIndex][LocalVectorIndex]*
                                                     FirstData[GlobalIndex][LocalVectorIndex];
                                         }
                                         size = Math.Sqrt(size);
                                         size = Math.Max(MinimumDistance, size);
                                         weight *= FirstScale/size;
                                     }

                                     var ValueFl = new double[PointVectorDimension];
                                     MatrixVectorProduct(ValueFl, Rotation, SecondData[GlobalIndex]);
                                     VectorSum(ValueFl, Translation, +1.0, ValueFl);
                                     VectorSum(ValueFl, ValueFl, -1.0, FirstData[GlobalIndex]);
                                     for (int LocalVectorIndex = 0;
                                          LocalVectorIndex < PointVectorDimension;
                                          LocalVectorIndex++)
                                         ValueFl[LocalVectorIndex] *= weight;


                                     for (int ipar = 0; ipar < Hotsun.npar; ipar++)
                                     {
                                         if (ipar < PointVectorDimension)
                                         {
                                             // Translation Derivatives
                                             for (int LocalVectorIndex = 0;
                                                  LocalVectorIndex < PointVectorDimension;
                                                  LocalVectorIndex++)
                                             {
                                                 DerivativeFl[LocalVectorIndex][ipar] = weight*
                                                                                        DerivTranslation[ipar][
                                                                                            LocalVectorIndex];
                                             }
                                         }
                                         else
                                         {
                                             var tempvector = new double[PointVectorDimension];
                                             MatrixVectorProduct(tempvector, DerivRotationScale[ipar],
                                                                 SecondData[GlobalIndex]);
                                             for (int LocalVectorIndex = 0;
                                                  LocalVectorIndex < PointVectorDimension;
                                                  LocalVectorIndex++)
                                             {
                                                 DerivativeFl[LocalVectorIndex][ipar] = weight*
                                                                                        tempvector[LocalVectorIndex];
                                             }
                                         }
                                     }
                                     for (int LocalVectorIndex = 0;
                                          LocalVectorIndex < PointVectorDimension;
                                          LocalVectorIndex++)
                                     {
                                         GenericManxcat.Accum(ThreadNo, ValueFl[LocalVectorIndex],
                                                              DerivativeFl[LocalVectorIndex]);
                                         if (GlobalIndex <= 0)
                                         {
                                             string deriv = "";
                                             for (int ipar = 0; ipar < Hotsun.npar; ipar++)
                                                 deriv += DerivativeFl[LocalVectorIndex][ipar].ToString("E4") + " ";
                                             // SALSAUtility.SALSAPrint(0, "Val " +  ValueFl[LocalVectorIndex].ToString("f4") + " Drv " + deriv);
                                             GenericManxcat.AccumDebug(GlobalIndex);
                                         }
                                         if (GlobalIndex > (SALSAUtility.PointCount_Global - 2))
                                             GenericManxcat.AccumDebug(GlobalIndex);
                                     }
                                 } // End loop over points in this thread
                             }); // End loop initialing Point dependent quantities

            // Add up MPI and Thread Contributions
            GenericManxcat.AddupChisqContributions(Solution);

            return violat;
        }

        // End Calcfg

        //  Convert angles to range -PI to PI
        public static double ProperAngle(double Angle)
        {
            double fudge = (Angle + Math.PI)/(2.0*Math.PI);
            return Angle - Math.Floor(fudge)*2.0*Math.PI;
        }

        // End ProperAngle(double Angle)

        public static void SetRotationMatrix(double[,] RotationMatrix, int FixedAxis, double Theta, bool Differentiate)
        {
            int AnotherAxis1 = (FixedAxis + 1)%PointVectorDimension;
            int itmp = (FixedAxis + 2)%PointVectorDimension;
            int AnotherAxis2 = itmp;
            if (AnotherAxis1 > AnotherAxis2)
            {
                AnotherAxis2 = AnotherAxis1;
                AnotherAxis1 = itmp;
            }

            if (Differentiate)
            {
                RotationMatrix[AnotherAxis1, AnotherAxis1] = - Math.Sin(Theta);
                RotationMatrix[AnotherAxis2, AnotherAxis2] = RotationMatrix[AnotherAxis1, AnotherAxis1];
                RotationMatrix[AnotherAxis1, AnotherAxis2] = Math.Cos(Theta);
                RotationMatrix[AnotherAxis2, AnotherAxis1] = -RotationMatrix[AnotherAxis1, AnotherAxis2];
            }
            else
            {
                if (PointVectorDimension == 3)
                    RotationMatrix[FixedAxis, FixedAxis] = 1.0;
                RotationMatrix[AnotherAxis1, AnotherAxis1] = Math.Cos(Theta);
                RotationMatrix[AnotherAxis2, AnotherAxis2] = RotationMatrix[AnotherAxis1, AnotherAxis1];
                RotationMatrix[AnotherAxis1, AnotherAxis2] = Math.Sin(Theta);
                RotationMatrix[AnotherAxis2, AnotherAxis1] = -RotationMatrix[AnotherAxis1, AnotherAxis2];
            }
        }

        // End SetRotationMatrix(double[,] RotationMatrix, int FixedAxis, double Theta, bool Differentiate)

        public static void SetTotalRotation(double[,] RTotal, double Scaling, double[,] Cosmic, double[,] R3,
                                            double[,] R2, double[,] R1)
        {
            var Rtemp1 = new double[PointVectorDimension,PointVectorDimension];
            var Rtemp2 = new double[PointVectorDimension,PointVectorDimension];
            RotateMatricies(Rtemp1, R2, R1);
            RotateMatricies(Rtemp2, R3, Rtemp1);
            RotateMatricies(RTotal, Cosmic, Rtemp2);
            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
            {
                for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                {
                    RTotal[LocalVectorIndex1, LocalVectorIndex2] *= Scaling;
                }
            }
        }

        // End SetTotalRotation(double[,] RTotal, double Scaling, double[,] Cosmic, double[,] R3, double[,] R2, double[,] R1)

        public static void SetTotalRotation(double[,] RTotal, double Scaling, double[,] Cosmic, double[,] R1)
        {
            RotateMatricies(RTotal, Cosmic, R1);
            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
            {
                for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                {
                    RTotal[LocalVectorIndex1, LocalVectorIndex2] *= Scaling;
                }
            }
        }

        // End SetTotalRotation(double[,] RTotal, double Scaling, double[,] Cosmic, double[,] R1)

        public static void RotateMatricies(double[,] RTotal, double[,] R2, double[,] R1)
        {
            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
            {
                for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                {
                    double tmp = 0.0;
                    for (int LocalVectorIndex3 = 0; LocalVectorIndex3 < PointVectorDimension; LocalVectorIndex3++)
                    {
                        tmp = tmp + R2[LocalVectorIndex1, LocalVectorIndex3]*R1[LocalVectorIndex3, LocalVectorIndex2];
                    }
                    RTotal[LocalVectorIndex1, LocalVectorIndex2] = tmp;
                }
            }
        }

        // End RotateMatricies(double[,] RTotal, double[,] R2, double[,] R1)

        public static void MatrixVectorProduct(double[] VTotal, double[,] R2, double[] V)
        {
            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
            {
                double tmp = 0.0;
                for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                {
                    tmp = tmp + R2[LocalVectorIndex1, LocalVectorIndex2]*V[LocalVectorIndex2];
                }
                VTotal[LocalVectorIndex1] = tmp;
            }
        }

        // End MatrixVectorProduct(double[] VTotal, double[,] R2, double[] V)

        public static void VectorSum(double[] VTotal, double[] V2, double sign, double[] V1)
        {
            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
            {
                VTotal[LocalVectorIndex1] = V2[LocalVectorIndex1] + sign*V1[LocalVectorIndex1];
            }
        }

        // End VectorSum(double[] VTotal, double[] V2, double sign, double[] V1)

        public static void InitializeParameters(Desertwind Solution, int CountStartingPoints)
        {
            // Inversion
            InvertSolution = false;
            if ((CountStartingPoints/2)*2 != CountStartingPoints)
                InvertSolution = true;

            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
            {
                for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                {
                    CurrentCosmicScaling[LocalVectorIndex1, LocalVectorIndex2] = 0.0;
                }
                CurrentCosmicScaling[LocalVectorIndex1, LocalVectorIndex1] = 1.0;
            }
            if (InvertSolution)
            {
                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
                {
                    CurrentCosmicScaling[LocalVectorIndex1, LocalVectorIndex1] = -1.0;
                }
                if (PointVectorDimension == 2)
                    CurrentCosmicScaling[0, 0] = 1.0;
            }
            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
            {
                for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                {
                    SavedCosmicScaling[CountStartingPoints][LocalVectorIndex1, LocalVectorIndex2] =
                        CurrentCosmicScaling[LocalVectorIndex1, LocalVectorIndex2];
                }
            }

            //  Scaling
            double Scale = FirstScale/SecondScale;
            if (ScaleOption == 0)
                Solution.param[ScalePosition][0] = Scale;
            else
            {
                ScaleA1 = 0.5*(ScaleAlpha + 1.0/ScaleAlpha);
                ScaleA2 = 0.5*Math.Abs(ScaleAlpha - 1.0/ScaleAlpha);
                Solution.param[ScalePosition][0] = Math.Asin((ScaleA1 - 1.0)/ScaleA2);
                ScaleA1 *= Scale;
                ScaleA2 *= Scale;
            }

            var Randobject = new Random();
            var RandomAngles = new double[3];
            if (SALSAUtility.MPI_Rank == 0)
            {
                for (int irand = 0; irand < 3; irand++)
                    RandomAngles[irand] = Randobject.NextDouble();
            }
            if (SALSAUtility.MPI_Size > 1)
            {
                SALSAUtility.StartSubTimer(SALSAUtility.MPIBROADCASTTiming);
                SALSAUtility.MPI_communicator.Broadcast(ref RandomAngles, 0);
                SALSAUtility.StopSubTimer(SALSAUtility.MPIBROADCASTTiming);
            }


            //  Translation and Rotation
            for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
            {
                Solution.param[LocalVectorIndex][0] = FirstMean[LocalVectorIndex] - Scale*
                                                      CurrentCosmicScaling[LocalVectorIndex, LocalVectorIndex]*
                                                      SecondMean[LocalVectorIndex];
                if (PointVectorDimension == 3)
                {
                    Solution.param[LocalVectorIndex + PointVectorDimension][0] = 0.0;
                    if (CountStartingPoints > 1)
                    {
                        Solution.param[LocalVectorIndex + PointVectorDimension][0] = Math.PI*
                                                                                     RandomAngles[LocalVectorIndex];
                    }
                }
            }
            if (PointVectorDimension == 2)
            {
                Solution.param[2][0] = 0.0;
                if (CountStartingPoints > 1)
                {
                    Solution.param[2][0] = Math.PI*RandomAngles[0];
                }
            }
        }

        // End InitializeParameters

        public static void Sequel()
        {
            /* Sequel of RotateManxcatMDS work */

            // At the moment do nothing.
            return;
        }

        public static double CalculateScale(double InputParm)
        {
            if (ScaleOption == 0)
            {
                return InputParm;
            }
            return ScaleA1 - ScaleA2*Math.Sin(InputParm);
        }

        public static double DifferentiateScale(double InputParm)
        {
            if (ScaleOption == 0)
            {
                return 1.0;
            }
            return (- ScaleA2*Math.Cos(InputParm));
        }

        public static void WriteRotateMDS(string fname, int OutputOption, double[][] param, double[][] perr)
        {
            int LocalVectorDimension = param[0].GetLength(0);
            if (LocalVectorDimension != perr[0].GetLength(0))
                SALSAUtility.SALSAError("Inconsistent Dimensions Labels " + LocalVectorDimension.ToString() + " Perr " +
                                        perr[0].GetLength(0).ToString());

            ManxcatSection Configuration = ManxcatCentral.Configuration;

            // set current cosmic scaling
            for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
            {
                for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                {
                    CurrentCosmicScaling[LocalVectorIndex1, LocalVectorIndex2] =
                        SavedCosmicScaling[Hotsun.BestChisqLoop][LocalVectorIndex1, LocalVectorIndex2];
                }
            }

            // Write SALSA Properties File Header
            if (OutputOption == 0)
            {
                int filetype = 2;
                if ((SALSAUtility.NumberFixedPoints > 0) ||
                    (SALSAUtility.NumberVariedPoints < SALSAUtility.NumberOriginalPoints))
                    filetype = 4;
                RotationFileProperties.GroupName = "MDS Run:" + Configuration.RunNumber.ToString();
                RotationFileProperties.FileGenerationType = filetype;
                RotationFileProperties.OriginalPointStartIndex = 0;
                RotationFileProperties.LocalPointStartIndex = 0;

                //  Rotation parameters
                RotationFileProperties.NumberRotationParameters = Hotsun.npar;
                RotationFileProperties.RotationParameters = "";
                for (int CountRotPars = 0; CountRotPars < Hotsun.npar; CountRotPars++)
                {
                    RotationFileProperties.RotationParameters += " " + param[CountRotPars][0].ToString("E4") + ",";
                }
                for (int Cosmicdiagonals = 0; Cosmicdiagonals < PointVectorDimension; Cosmicdiagonals++)
                {
                    RotationFileProperties.RotationParameters += " " +
                                                                 CurrentCosmicScaling[Cosmicdiagonals, Cosmicdiagonals].
                                                                     ToString("E4") + ",";
                }
                // Comment should have key features of Run
                string OldComment = RotationFileProperties.Comment;
                if (OldComment != "")
                    OldComment += "\n";
                OldComment += "MDS ROTATION Run " + Configuration.RunNumber.ToString() + " Start Time " +
                              SALSAUtility.startTime.ToLocalTime() + " *** ";
                OldComment += "\n RotationOption " + Configuration.RotationOption + " First Scale " +
                              FirstScale.ToString("E4") + " Second Scale " + SecondScale.ToString()
                              + "\n First Mean";
                for (int PointIndex = 0; PointIndex < PointVectorDimension; PointIndex++)
                {
                    OldComment += " " + FirstMean[PointIndex].ToString("E4");
                }
                OldComment += " Second Mean";
                for (int PointIndex = 0; PointIndex < PointVectorDimension; PointIndex++)
                {
                    OldComment += " " + SecondMean[PointIndex].ToString("E4");
                }
                RotationFileProperties.Comment = OldComment;

                SALSA_Properties.WriteDataPointFile(fname, ManxcatCentral.Configuration.Write2Das3D, "all ",
                                                    RotationFileProperties,
                                                    RotationPointProperties, SALSAUtility.NumberOriginalPoints);
                return;
            }

            //  Option OutputOption = 1 -- write ROTATED Second Data
            try
            {
                StreamWriter sw = null;
                if (!string.IsNullOrEmpty(fname))
                {
                    sw = new StreamWriter(fname, false, Encoding.UTF8);
                }
                if (sw != null)
                {
                    double[] FullTranslation;
                    double[,] FullRotation;
                    SetupFinalTransforamation(param, out FullTranslation, out FullRotation);

                    for (int UsedDataPoint = 0; UsedDataPoint < SALSAUtility.PointCount_Global; UsedDataPoint++)
                    {
                        var Vector = new double[PointVectorDimension];
                        MatrixVectorProduct(Vector, FullRotation, SecondData[UsedDataPoint]);
                        VectorSum(Vector, FullTranslation, +1.0, Vector);

                        string Coordinates = "";
                        int SingleCluster = 1;

                        for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
                        {
                            Coordinates += Vector[LocalVectorIndex].ToString("E4") + "\t";
                        }
                        sw.WriteLine(
                            String.Format((UsedDataPoint + 1).ToString() + "\t" + Coordinates + SingleCluster.ToString()));
                    }
                }

                sw.Flush();
                sw.Close();
            }
            catch (Exception e)
            {
                Console.WriteLine("Failed writing data" + e);
                throw (e);
            }
        }

        // End WriteRotateMDS

        // Scale is now sqrt of average squared distance from mean of distribution
        public static void ProcessPointData(SALSADataPointProperties[] DataFile, double[][] NewPointArray,
                                            double[] PointMean, out double Scale)
        {
            for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
            {
                PointMean[LocalVectorIndex] = 0.0;
            }
            for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
            {
                int OriginalIndex = SALSAUtility.UsedPointtoOriginalPointMap[GlobalPointIndex];
                PointMean[0] += DataFile[OriginalIndex].x;
                PointMean[1] += DataFile[OriginalIndex].y;

                NewPointArray[GlobalPointIndex][0] = DataFile[OriginalIndex].x;
                NewPointArray[GlobalPointIndex][1] = DataFile[OriginalIndex].y;
                if (PointVectorDimension > 2)
                {
                    PointMean[2] += DataFile[OriginalIndex].z;
                    NewPointArray[GlobalPointIndex][2] = DataFile[OriginalIndex].z;
                }
            }
            for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
            {
                PointMean[LocalVectorIndex] /= SALSAUtility.PointCount_Global;
            }

            Scale = 0.0;
            for (int GlobalPointIndex = 0; GlobalPointIndex < SALSAUtility.PointCount_Global; GlobalPointIndex++)
            {
                double tmp1 = 0.0;
                for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
                {
                    double tmp2 = NewPointArray[GlobalPointIndex][LocalVectorIndex] - PointMean[LocalVectorIndex];
                    tmp1 += tmp2*tmp2;
                }
                Scale += tmp1;
            }
            Scale = Math.Sqrt(Scale/SALSAUtility.PointCount_Global);
        }

        // End ProcessPointData

        public static void SetupFinalTransforamation(double[][] param, out double[] Translation, out double[,] Rotation)
        {
            // Calculate Transformation without derivations

            double Scale = CalculateScale(param[ScalePosition][0]);

            //  Each Value is weight * (Scale*CosmicScaling*Rotation*Second Vector + Translation - First Vector)
            //  Each Derivative is weight * ( Rotation Term*Second Vectot + Translation Term )
            //  Rotation, Rotation Term, Translation, Translation Term are independent of Point
            //  weight depens on Option Used 

            Translation = new double[PointVectorDimension];
            Rotation = new double[PointVectorDimension,PointVectorDimension];
            var SubRotation = new double[NumberSubRotations][,];

            for (int LocalVectorIndex = 0; LocalVectorIndex < PointVectorDimension; LocalVectorIndex++)
            {
                Translation[LocalVectorIndex] = param[LocalVectorIndex][0]; // Translations
            }

            for (int LocalRotationIndex = 0; LocalRotationIndex < NumberSubRotations; LocalRotationIndex++)
            {
                SubRotation[LocalRotationIndex] = new double[PointVectorDimension,PointVectorDimension];
                for (int LocalVectorIndex1 = 0; LocalVectorIndex1 < PointVectorDimension; LocalVectorIndex1++)
                {
                    for (int LocalVectorIndex2 = 0; LocalVectorIndex2 < PointVectorDimension; LocalVectorIndex2++)
                    {
                        SubRotation[LocalRotationIndex][LocalVectorIndex1, LocalVectorIndex2] = 0.0;
                    }
                }
                SetRotationMatrix(SubRotation[LocalRotationIndex], LocalRotationIndex,
                                  param[PointVectorDimension + LocalRotationIndex][0], false);
            }
            SetTotalRotation(Rotation, Scale, CurrentCosmicScaling, SubRotation[2], SubRotation[1], SubRotation[0]);
        }

        // End SetupFinalTransformation
    }
}