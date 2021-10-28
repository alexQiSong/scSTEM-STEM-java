package edu.cmu.cs.sb.stem;

import edu.cmu.cs.sb.chromviewer.*;
import edu.cmu.cs.sb.core.*;
import java.util.*;
import javax.swing.*;
import java.awt.event.*;

/**
 * Class encapsulates information on the comparison of two STEM data sets
 * focusing on identifying pairs of profiles with significant intersections 
 * between them.
 */
public class CompareInfo 
{
    /**
     * The main reference data set
     */
    STEM_DataSet origset;

    /**
     * The comparison data set
     */
    STEM_DataSet comparesetfmnel;

    /**
     * The p-value threshold for a significant interestection
     */
    double dmaxpval;

    /**
     * The mininum number of genes in a significant intersection
     */
    int nminnumgenes;

    /**
     * Arrays of CompareInfoRow about significant intersections for one priofile
     */
    ArrayList sigrows;

    /**
     * The version of sigrows transposed when swapping which data set is positioned
     * as the reference and which is the comparison
     */
    ArrayList sigrowsswap;
    
    
    ArrayList sigClusterPairs;

    /**
     * The maingui window for the comparison data
     */
    MAINGUI2 compareframe;

    /**
     * true if 2 datasets have the same number of time points, false otherwise
     */
    boolean ColFlag;

    /**
     * Use clusters or model profiles for comparison
     */
    boolean ClusterFlag;

    static int[] clusterID1;
    static int[] clusterID2;

    /**
     * Launches a maingui for the comparison data and finds the significant pairs of intersections
     */
    public CompareInfo(STEM_DataSet origset, String szcompare, 
		       String szmaxpval, String szminnumgenes,
                       Vector repeatnames, boolean balltime, ChromFrame cf, boolean ClusterFlag) throws Exception
    {

       this.origset = origset;
       this.ClusterFlag = ClusterFlag;


       dmaxpval = Double.parseDouble(szmaxpval);
       nminnumgenes = Integer.parseInt(szminnumgenes);
       sigrows = new ArrayList();
       sigClusterPairs = new ArrayList();
       //signames = new ArrayList();

       //builds a STEM comparison data set for this new data set, using the same parameters as the original set
       comparesetfmnel = ST.buildset(
		     cf.genomeParser.szchromval,
		    origset.tga.szxrefval,
                    szcompare, origset.tga.szGoFile, 
                    origset.tga.szGoCategoryFile, origset.nmaxmissing,
                    origset.dthresholdvalue, origset.dmincorrelation,
                    origset.dminclustdist, origset.alpha,
                    origset.dpercentileclust, origset.nmaxchange,
                    origset.nmaxprofiles, origset.dmaxcorrmodel,
                    origset.nsamplesgene, origset.tga.nsamplespval, 
                    origset.nsamplesmodel, origset.tga.nmingo, 
                    origset.nfdr, origset.tga.nmingolevel,
                    origset.tga.szextraval, balltime, repeatnames,
		    origset.btakelog,origset.bspotincluded,
                    origset.badd0,origset.tga.szcategoryIDval,
                    origset.tga.szevidenceval,origset.tga.sztaxonval,
                    origset.tga.bpontoval,origset.tga.bcontoval,origset.tga.bfontoval,
                    origset.tga.brandomgoval,origset.bkmeans,
                    origset.bmaxminval,origset.ballpermuteval,
		    origset.tga.szorganismsourceval, origset.tga.szxrefsourceval,
		    cf.genomeParser.szchromsourceval);

      
       if (comparesetfmnel.numcols != origset.numcols)
       {
           this.ColFlag = false;
       }
       else{
           this.ColFlag = true;
       }
 
 
       //launches an interface for this comparison dataset
       compareframe = new MAINGUI2(comparesetfmnel);
       compareframe.cf = cf;
       edu.umd.cs.piccolo.PCanvas.CURRENT_ZCANVAS = null;
 
       compareframe.setLocation(25,60);
       compareframe.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
       compareframe.addWindowListener(new WindowAdapter() {
          public void windowClosing(WindowEvent we) {
		    compareframe.closeSortWindows();
		}
                 public void windowClosed(WindowEvent we) {
                     Thread t = new Thread (new Runnable() { 
                                              public void run() { 
		                                 System.gc();  
                                              } 
                                          } 
                                    ); 
                     t.start(); 
		 }
	      });
	 
        //compareframe.setVisible(true);
        compareframe.setVisible(false);
        

        compareframe.dispose();
        
        
        if (!ClusterFlag)
        {
            findsigcompare();
        }
        else
        {
            FindSigClusterCompare();
        }

    }  

    /**
     * Returns the number of common elements in hm1 and hm2
     */
    private int countintersect(HashMap hm1, HashMap hm2)
    {
	Iterator em1 = hm1.keySet().iterator();
        int ncount = hm2.size();
        while (em1.hasNext())
	{
	    String sz1 = (String) em1.next();
            if (hm2.get(sz1) == null)
	    {
		ncount++;
	    }
	}
        return ncount;
    }

    /**
     * Finds significant pairs of profiles from the two data sets interms of an intersection
     * exceeding nminnumgenes and p-value less than dmaxpval
     */
    public void findsigcompare()
    {

        int nsize1, nsize2;
        ArrayList[] sigprofilesswapArray = new ArrayList[comparesetfmnel.modelprofiles.length];
        ArrayList[] sigprofilesswapnames = new ArrayList[comparesetfmnel.modelprofiles.length];
        double[] mincorrswap = new double[comparesetfmnel.modelprofiles.length];
        double[] minpvalswap = new double[comparesetfmnel.modelprofiles.length];

        for (int nindex = 0; nindex < sigprofilesswapArray.length; nindex++)
        {
	   sigprofilesswapArray[nindex] = new ArrayList();
	   sigprofilesswapnames[nindex] = new ArrayList();
           mincorrswap[nindex] = 1;
           minpvalswap[nindex] = 1;
	}

        int numintersect = countintersect(origset.tga.htGeneNames, comparesetfmnel.tga.htGeneNames);
	System.out.println("Number of genes in the union of the original and comparison set is "+numintersect);

        for (int nprofileindex = 0; nprofileindex < origset.modelprofiles.length; nprofileindex++)
	{
            HashSet htprofilegenes = new HashSet();
            nsize1 =origset.profilesAssigned[nprofileindex].size();
            for (int ngeneindex = 0; ngeneindex < nsize1; ngeneindex++)
	    {
                int nrealgeneindex = 
                          ((Integer) origset.profilesAssigned[nprofileindex].get(ngeneindex)).intValue();
                htprofilegenes.add(origset.genenames[nrealgeneindex]);
	    }
               
            ArrayList sigprofiles = new ArrayList();
            ArrayList sigprofilesnames = new ArrayList();   
            double dminpval = 1;
            double dmincorr = 1;
            for (int nprofileindex2 = 0; nprofileindex2 < comparesetfmnel.modelprofiles.length; nprofileindex2++)
	    {
               HashSet intersectnames = new HashSet();  
               double dmatch = 0;
               nsize2 = comparesetfmnel.profilesAssigned[nprofileindex2].size();

	       for (int ngeneindex = 0; ngeneindex < nsize2; ngeneindex++)
	       {
                  int nrealgeneindex = 
                        ((Integer) comparesetfmnel.profilesAssigned[nprofileindex2].get(ngeneindex)).intValue();
                  String szgene = comparesetfmnel.genenames[nrealgeneindex];
                  if (htprofilegenes.contains(szgene))
		  {
                      intersectnames.add(szgene);
                      dmatch +=1.0/comparesetfmnel.bestassignments[nrealgeneindex].size();
                  }
	       }

              double dpval =StatUtil.hypergeometrictail((int) Math.ceil(dmatch-1),nsize1,
			   (int)Math.floor(numintersect-nsize1),nsize2);
               if ((dmatch >=nminnumgenes)&&(dpval <= dmaxpval))
	       {
		  if (dpval < dminpval)
		  {
		     dminpval = dpval;
		  }

                  if (dpval < minpvalswap[nprofileindex2])
		  {
		      minpvalswap[nprofileindex2] = dpval;
		  }
                 
                 double dcorrval = 1;
                 if (this.ColFlag == true) 
                 {
                      dcorrval = Util.correlation(origset.modelprofiles[nprofileindex],
						      comparesetfmnel.modelprofiles[nprofileindex2]);
		      if (dcorrval < dmincorr)
		      {
                          dmincorr = dcorrval; 
		      }
		      
                     if (dcorrval < mincorrswap[nprofileindex2])
		      {
		           mincorrswap[nprofileindex2] = dcorrval;
		      }	      
                 }
                 else{
                 
		      dcorrval = 1;
                 }




		  
		  
                  sigprofiles.add(new CompareInfoRec(dmatch, dpval, dcorrval, nprofileindex2, intersectnames));
                  sigprofilesswapArray[nprofileindex2].add(new CompareInfoRec(dmatch,dpval,dcorrval,
									  nprofileindex,intersectnames));
	       }
	    }

            if (sigprofiles.size()>=1)
	    {
                sigrows.add(new CompareInfoRow(nprofileindex, dmincorr, dminpval,sigprofiles));
	    }
	}

        sigrowsswap = new ArrayList();
        for (int nprofileindex2 = 0; nprofileindex2 < sigprofilesswapArray.length; nprofileindex2++)
	{
            if  (sigprofilesswapArray[nprofileindex2].size() >= 1)
	    {
		sigrowsswap.add(new CompareInfoRow(nprofileindex2,mincorrswap[nprofileindex2],
			      minpvalswap[nprofileindex2],sigprofilesswapArray[nprofileindex2]));
	    }
	}
    }


    /**
     * Information on an intersection match: the number in the match, the p-value,
     * the correlation of the profiles, the profile ID, and the genes in the intersection.
     */
    static class CompareInfoRec
    {
	double dmatch;
        double dpval;
        double dcorrval;
        int nprofile;
        HashSet inames;

        CompareInfoRec(double dmatch, double dpval, double dcorrval, int nprofile, HashSet inames)
	{
	    this.dmatch = dmatch;
            this.dpval  = dpval;
            this.dcorrval = dcorrval;
            this.nprofile = nprofile;
            this.inames  = inames;
	}
    }

    /**
     *For a row of intersections: the profile ID from one data set, the list of significant intersections in
     * signprofiles, and the minimum correlation and p-value with any of these significant intersections
     */
   static class CompareInfoRow
   {
       int nprofile;
       double dmincorr;
       double dminpval;
       ArrayList sigprofiles;

       CompareInfoRow(int nprofile,double dmincorr, double dminpval, ArrayList sigprofiles)
       {
	   this.dmincorr = dmincorr;
	   this.dminpval = dminpval;
           this.nprofile = nprofile;
           this.sigprofiles = sigprofiles;
       }
   }
   


    /**
     * Information on an intersection match of 2 clusters: the number in the match, the p-value,
     * the correlation of the profiles, the profile ID, and the genes in the intersection.
     */
    public static final class sigClusterPair
    {
	public double[] dmatch;
        public double[] dpvals;
        public double dcorrval;
        public ArrayList cluster1;
        public ArrayList rightClusters;
        public ArrayList inames;

        sigClusterPair(double[] dmatch, double[] dpvals, double dcorrval, ArrayList cluster1, ArrayList rightClusters, ArrayList intersectNameList)
	{
	    this.dmatch = dmatch;
            this.dpvals  = dpvals;
            this.dcorrval = dcorrval;
            this.cluster1 = cluster1;   //ArrayList
            this.rightClusters = rightClusters; //Arraylist of clusters
            this.inames  = intersectNameList; // ArrayList of intersectnames(HashSet)
	}
    }


   /**
     * Finds significant pairs of CLUSTERS from the two data sets interms of an intersection
     * exceeding nminnumgenes and p-value less than dmaxpval
     */
    public void FindSigClusterCompare()
    {
        System.out.println("Comparing clusters...");

        this.clusterID1 = new int[origset.modelprofiles.length];
        this.clusterID2 = new int[comparesetfmnel.modelprofiles.length];
        
        
        //get profile cluster IDs
        int numclusters2 = comparesetfmnel.clustersofprofilesnum.size();
        for (int nindex = 0; nindex < comparesetfmnel.modelprofiles.length; nindex++)
	{
           this.clusterID2[nindex] = numclusters2;  
	}     

        for (int ncluster= 0; ncluster< comparesetfmnel.clustersofprofilesnum.size(); ncluster++)
        {
	   ArrayList profilesInCluster = (ArrayList) comparesetfmnel.clustersofprofilesnum.get(ncluster);
           for (int nprofileindex = 0; nprofileindex < profilesInCluster.size(); nprofileindex++)
	   {
	      STEM_DataSet.ProfileRec theProfile = (STEM_DataSet.ProfileRec) profilesInCluster.get(nprofileindex);
              this.clusterID2[theProfile.nprofileindex] = ncluster;           
	   }
	}  
        
        int numclusters1 = origset.clustersofprofilesnum.size();
        for (int nindex = 0; nindex < origset.modelprofiles.length; nindex++)
	{
           this.clusterID1[nindex] = numclusters1;
	}     

        for (int ncluster= 0; ncluster< origset.clustersofprofilesnum.size(); ncluster++)
        {
	   ArrayList profilesInCluster = (ArrayList) origset.clustersofprofilesnum.get(ncluster);
           for (int nprofileindex = 0; nprofileindex < profilesInCluster.size(); nprofileindex++)
	   {
	      STEM_DataSet.ProfileRec theProfile = (STEM_DataSet.ProfileRec) profilesInCluster.get(nprofileindex);
              this.clusterID1[theProfile.nprofileindex] = ncluster;             
	   }
	}  
        // end of getting clusters
        

        
        // get max cluster IDs, which indicates profiles not in clusters
        int maxID1 = 0;
        for (int i=0; i<this.clusterID1.length;i++)
        {
            if (this.clusterID1[i]>maxID1)
            {
                maxID1 = this.clusterID1[i];
            }
        }
        int maxID2 = 0;
        for (int i=0; i<this.clusterID2.length;i++)
        {
            if (this.clusterID2[i]>maxID2)
            {
                maxID2 = this.clusterID2[i];
            }
        }   


      
      
        ArrayList[] SigClustersSwapArray = new ArrayList[numclusters2];
        ArrayList[] SigClustersSwapNames = new ArrayList[numclusters2];
    
        
        double[] mincorrswap = new double[numclusters2];
        double[] minpvalswap = new double[numclusters2];


        for (int nindex = 0; nindex < numclusters2; nindex++)
        {
	     SigClustersSwapArray[nindex] = new ArrayList();
	     SigClustersSwapNames[nindex] = new ArrayList();
             mincorrswap[nindex] = 1;
             minpvalswap[nindex] = 1;
	}



        int numintersect = countintersect(origset.tga.htGeneNames, comparesetfmnel.tga.htGeneNames);
	System.out.println("Number of genes in the union of the original and comparison set is "+numintersect);


        for (int i = 0; i < origset.clustersofprofilesnum.size(); i++) //loop over orignal set clusters
	{
	    ArrayList pvals = new ArrayList(); 
            int nsize1=0;
            int nsize2=0;
            HashSet htprofilegenes = new HashSet();
            ArrayList cluster = (ArrayList) origset.clustersofprofilesnum.get(i);
            ArrayList intersectNameList = new ArrayList(); // a list of intersect names in a row.
            ArrayList comparedCluster1 = (ArrayList)(origset.clustersofprofilesnum).get(i);
            
            for (int p=0; p<cluster.size(); p++ ) // for each profile, get genes
            {
                STEM_DataSet.ProfileRec profile = (STEM_DataSet.ProfileRec)((STEM_DataSet.ProfileRec)cluster.get(p)).clone();
                int nprofileindex = profile.nprofileindex;
                int psize1 =origset.profilesAssigned[nprofileindex].size();
                nsize1 += psize1;
                for (int ngeneindex = 0; ngeneindex < psize1; ngeneindex++)
	        {
                    int nrealgeneindex = 
                              ((Integer) origset.profilesAssigned[nprofileindex].get(ngeneindex)).intValue();
                    htprofilegenes.add(origset.genenames[nrealgeneindex]);
	        }    
            }
            
            ArrayList rightClusters = new ArrayList(); // list for storing clusters and profiles compared (RHS. of the UI)
            
            double dminpval = 1;
            double dmincorr = 1;
            ArrayList dmatchlist = new ArrayList();
            //double dmatch = 0;    
            for (int j = 0; j < comparesetfmnel.clustersofprofilesnum.size(); j++) // loop over comparison set clusters 
	    {
                double dmatch = 0;
                nsize2=0;
                HashSet intersectnames = new HashSet();  
                ArrayList cluster2 = (ArrayList) comparesetfmnel.clustersofprofilesnum.get(j);
                for (int p2=0; p2<cluster2.size(); p2++)
                {
                    STEM_DataSet.ProfileRec profile2 = (STEM_DataSet.ProfileRec)((STEM_DataSet.ProfileRec)cluster2.get(p2)).clone();
                    int nprofileindex2 = profile2.nprofileindex;
                    int psize2 =comparesetfmnel.profilesAssigned[nprofileindex2].size();
                    nsize2 += psize2;
                    for (int ngeneindex2 = 0; ngeneindex2 < psize2; ngeneindex2++)
	            {
                        int nrealgeneindex2 = 
                              ((Integer) comparesetfmnel.profilesAssigned[nprofileindex2].get(ngeneindex2)).intValue();
                        String szgene = comparesetfmnel.genenames[nrealgeneindex2];
                        if (htprofilegenes.contains(szgene))
                        {
                            intersectnames.add(szgene);
                            dmatch +=1.0/comparesetfmnel.bestassignments[nrealgeneindex2].size();
                            //dmatch += 1.0;
                        }
                        
	            }    
                }



                            
                double dpval =StatUtil.hypergeometrictail((int) Math.ceil(dmatch-1),nsize1,
			   (int)Math.floor(numintersect-nsize1),nsize2);
                if ((dmatch >=nminnumgenes)&&(dpval <= dmaxpval))
  	        {
	            if (dpval < dminpval)
		    {
		        dminpval = dpval;
		    }

                    if (dpval < minpvalswap[j])
		    {
		        minpvalswap[j] = dpval;
		    }
                 
                   double dcorrval = 1; // not include correlations now

                   
                   ArrayList comparedCluster2 = (ArrayList)(comparesetfmnel.clustersofprofilesnum).get(j);
                   rightClusters.add(comparedCluster2);
                   for (int p2=0; p2<cluster2.size(); p2++)  // add intersect names for each profile in the cluster
                   {
                       intersectNameList.add(intersectnames);
                       pvals.add(dpval);  // add pvals for each profile in the cluster
                       dmatchlist.add(dmatch); // add dmatch for each profile
                   }
                      
	        }
            } // end of comparison cluster comparison




            ArrayList sigprofiles = new ArrayList();
	   
            double dminpvalp = 1;
            double dmincorrp = 1;
	    
            for (int nprofileindex2 = 0; nprofileindex2 < comparesetfmnel.modelprofiles.length; nprofileindex2++)
            {   
            // compare the cluster with profiles
                if (this.clusterID2[nprofileindex2] != maxID2)
                {
                    continue; //its a profile in a cluster, skip since already compared.
                }
                HashSet intersectnamesp = new HashSet();  
                double dmatchp = 0;
                nsize2 = comparesetfmnel.profilesAssigned[nprofileindex2].size();
                    
                //intersectNameList.add(intersectnamesp);
                
                for (int ngeneindex = 0; ngeneindex < nsize2; ngeneindex++)
	        {
                    int nrealgeneindex = 
                        ((Integer) comparesetfmnel.profilesAssigned[nprofileindex2].get(ngeneindex)).intValue();
                    String szgene = comparesetfmnel.genenames[nrealgeneindex];
                    if (htprofilegenes.contains(szgene))
                    {
                        intersectnamesp.add(szgene);
                        dmatchp +=1.0/comparesetfmnel.bestassignments[nrealgeneindex].size();
                        //dmatchp += 1.0;
                    }
                }
                
                
                double dpval =StatUtil.hypergeometrictail((int) Math.ceil(dmatchp-1),nsize1,
			   (int)Math.floor(numintersect-nsize1),nsize2);
                if ((dmatchp >=nminnumgenes)&&(dpval <= dmaxpval))
                {   // if significant 
                    if (dpval < dminpval)
                    {
                        dminpval = dpval;
                    }
                    dmatchlist.add(dmatchp);
                    pvals.add(dpval);
                    intersectNameList.add(intersectnamesp);
                    ArrayList profileCluster = new ArrayList();
                    STEM_DataSet.ProfileRec dummyProfile = new STEM_DataSet.ProfileRec(nprofileindex2, 0,0,dpval);
                    profileCluster.add(dummyProfile);
                    rightClusters.add(profileCluster);  // make this significant profile a cluster and add to the RHS. clusters list
                    
                } 
            } //end of cluster vs profile loop
            
            if (comparedCluster1.size()!=0 && rightClusters.size()!=0)
            {
                double[] dpvals = new double[pvals.size()];
                for (int index=0; index<pvals.size();index++)
                {
                    dpvals[index]= (double) pvals.get(index);
                }
                
                double[] dmatchdouble = new double[dmatchlist.size()];
                for (int matchindex=0; matchindex<dmatchlist.size();matchindex++)
                {
                    dmatchdouble[matchindex] = (double) dmatchlist.get(matchindex);
                }

                
                sigClusterPair pair = new sigClusterPair(dmatchdouble, dpvals,1,comparedCluster1, rightClusters, intersectNameList);
                sigClusterPairs.add(pair);
            }
            
        
        } // end of outer loop

        
    }






}
