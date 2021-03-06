package edu.cmu.cs.sb.stem;

import edu.cmu.cs.sb.chromviewer.*;
import edu.cmu.cs.sb.core.*;
import edu.umd.cs.piccolo.PCanvas;
import edu.umd.cs.piccolo.PNode;
import edu.umd.cs.piccolo.event.PZoomEventHandler;
import edu.umd.cs.piccolo.event.PBasicInputEventHandler;
import edu.umd.cs.piccolo.event.PInputEvent;
import edu.umd.cs.piccolo.nodes.PPath;
import edu.umd.cs.piccolox.PFrame;
import edu.umd.cs.piccolo.nodes.PText;
import edu.umd.cs.piccolo.nodes.PImage;
import java.util.*;
import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import java.net.*;


/**
 * Class for the window that compares profiles
 */
public class CompareGui extends PFrame
{
    /**
     * Screen width 
     */
    static int SCREENWIDTH = 800;

    /**
     * Screen height
     */
    static int SCREENHEIGHT = 600;

    /**
     * Space between edges and interface
     */
    static int BUFFER = 135;

    /**
     * Offset for the caption parameter
     */
    static int CAPOFFSET = 20;

    /**
     * constant for sorting profile pair by ID
     */
    static int SORTID = 0;

    /**
     * constant for sorting profile pair by significance intersection
     */
    static int SORTSIG = 1;

    /**
     * constant for sorting profile pair by suprise
     */
    static int SORTSUPRISE = 2;

    /**
     * The canvas on which the signifcant profile pairs are displayed
     */
    private PCanvas canvas;

    /**
     * The CompareInfo that this gui display 
     */
    CompareInfo theCompareInfo = null;

    /**
     * Array of the profiles associated with the original dataset
     */
    private PNode[] profilenodes;

    /**
     * specifies how to sort profile pairs
     */ 
    private int nsort= 0;

    /**
     * If the original and comparison data sets should be swapped from their original positions
     */
    private boolean bswap = false;

    /**
     * Used to synchronize swapping actions
     */
    private Object swaplock = new Object();

    /**
     * Label for the other data set in the layout whether it is the comparison or the original
     */
    private String szother;

    /**
     * String for labels indicating if we have profile pairs or captions
     */
    private String szprofileclusterCAP;

    /**
     * thegeneplotpanel associated with the original data set
     */
    private GenePlotPanel thegeneplotpanel;

    /**
     * Interface options corresponding to data set indexed as rows in the comparison
     */
    private GenePlotPanel rowplotpanel;


    private GenePlotPanel leftplotpanel;
    

    /**
     * Interface options corresponding to data set currently indexed as columns in the comparison
     */
    private GenePlotPanel colplotpanel;


    private GenePlotPanel rightplotpanel;

    /**
     * The chromosome viewer object
     */
    private ChromFrame cf;

    /**
     * Use clusters or model profiles for comparison
     */
    public boolean ClusterFlag;


    /**
     * Creates a comparison intersection interface backed on theCompareInfo
     */
    public CompareGui(CompareInfo theCompareInfo,PNode[] profilenodes,
                      GenePlotPanel thegeneplotpanel,ChromFrame cf,boolean ClusterFlag) throws Exception
    {
       super("Comparison - Significant Intersections",false,null);

       this.theCompareInfo = theCompareInfo;
       this.profilenodes = profilenodes;
       this.thegeneplotpanel = thegeneplotpanel;
       this.cf = cf;
       this.ClusterFlag = ClusterFlag;
       
       //this.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);//exit the app when close this window
       
       if (theCompareInfo.origset.bkmeans)
       {
          szprofileclusterCAP = "Clusters";
       }
       else
       {
	  szprofileclusterCAP = "Profiles";
       }
    }

    /**
     * Sorts rows based on reference profile most significant intersection p-value, then profile ID
     */
    public static class SigRowComparator implements Comparator
    {
	public int compare(Object o1, Object o2)
	{
            CompareInfo.CompareInfoRow cro1 = (CompareInfo.CompareInfoRow) o1;
            CompareInfo.CompareInfoRow cro2 = (CompareInfo.CompareInfoRow) o2;

            if (cro1.dminpval < cro2.dminpval)
                return -1;
            else if (cro1.dminpval > cro2.dminpval)
		return 1;
            else if (cro1.nprofile < cro2.nprofile)
                return -1;
            else if (cro1.nprofile > cro2.nprofile)
                return 1;

            return 0;
	}
    }

    /**
     * Sorts individual profiles by p-value with reference profile in row, then profile ID 
     */
    public static class SigComparator implements Comparator
    {
	public int compare(Object o1, Object o2)
	{
	    CompareInfo.CompareInfoRec co1 = (CompareInfo.CompareInfoRec) o1;
	    CompareInfo.CompareInfoRec co2 = (CompareInfo.CompareInfoRec) o2;
            if (co1.dpval < co2.dpval)
                return -1;
            else if (co1.dpval > co2.dpval)
		return 1;
            else if (co1.nprofile < co2.nprofile)
                return -1;
            else if (co1.nprofile > co2.nprofile)
                return 1;

            return 0;
	}
    }

    /**
     * Sorts rows based on minimum correlation with reference profile, then profile ID
     */
    public class SupriseComparator implements Comparator
    {
        int norigprofile;
        SupriseComparator(int norigprofile)
	{
	    this.norigprofile = norigprofile;
	}

	public int compare(Object o1, Object o2)
	{
	    CompareInfo.CompareInfoRec co1 = (CompareInfo.CompareInfoRec) o1;
	    CompareInfo.CompareInfoRec co2 = (CompareInfo.CompareInfoRec) o2;
            double dcorr1 = Util.correlation(
                                theCompareInfo.origset.modelprofiles[norigprofile],
				theCompareInfo.comparesetfmnel.modelprofiles[co1.nprofile]);

            double dcorr2 = Util.correlation(
                                theCompareInfo.origset.modelprofiles[norigprofile],
				theCompareInfo.comparesetfmnel.modelprofiles[co2.nprofile]);

            if (dcorr1 < dcorr2)
                return -1;
            else if (dcorr1 > dcorr2)
		return 1;
            else if (co1.nprofile < co2.nprofile)
                return -1;
            else if (co1.nprofile > co2.nprofile)
                return 1;

            return 0;
	}
    }


    /**
     * Sorts rows based on minimum correlation with reference profile, then profile ID
     */
    public static class SupriseRowComparator implements Comparator
    {
	public int compare(Object o1, Object o2)
	{
            CompareInfo.CompareInfoRow cro1 = (CompareInfo.CompareInfoRow) o1;
            CompareInfo.CompareInfoRow cro2 = (CompareInfo.CompareInfoRow) o2;

            if (cro1.dmincorr < cro2.dmincorr)
                return -1;
            else if (cro1.dmincorr > cro2.dmincorr)
		return 1;
            else if (cro1.nprofile < cro2.nprofile)
                return -1;
            else if (cro1.nprofile > cro2.nprofile)
                return 1;

            return 0;
	}
    }
    
    /**
     * sets the size of the PFrame
     */
    public void beforeInitialize()
    {
       setSize(SCREENWIDTH,SCREENHEIGHT);
    }

    /**
     * Simply calls drawcomparemain
     */
    public void initialize()
    {
	drawcomparemain();
    }

    /**
     * Lays out the comparison panel
     */
    public void drawcomparemain()
    {	    
        if (ClusterFlag) //compare clusters condition
        {
            System.out.println("Draw compare clusters UIs.");
            canvas = getCanvas();
            canvas.getLayer().removeAllChildren();
            int nmaxcol = 0;
            float frecwidth = (float) 150.0;
            int OFFSET = 0;
            ArrayList sigClusterPairs;
            final STEM_DataSet leftset, rightset;
            PNode[] leftnodes, rightnodes;
            String sztextleft, sztextright;


            sigClusterPairs = theCompareInfo.sigClusterPairs;
            leftset = theCompareInfo.origset;
            rightset = theCompareInfo.comparesetfmnel;
            leftnodes = profilenodes;
            rightnodes = theCompareInfo.compareframe.profilenodes;
            rightplotpanel = theCompareInfo.compareframe.thegeneplotpanel;
            leftplotpanel = thegeneplotpanel;
            sztextleft = ST.szDataFileDEF.split("_")[0];
            sztextright = ST.szCompareDEF.split("_")[0];

            szother = ST.szDataFileDEF.split("_")[0];

            PText comparetext = new PText(sztextright);
            comparetext.setFont(new Font("times",Font.PLAIN,14));
            double textXOffset1 = 0;

            PText originaltext = new PText(sztextleft);
            originaltext.setFont(new Font("times",Font.PLAIN,14));
            double textXOffset2 = 0;
            double scale1 = 0;
            double scale2 = 0;
            
            int nsigpairs = sigClusterPairs.size();
            
            
            ArrayList rhss = new ArrayList(); // rhs profiles of all rows
            
            
            
            for (int nrow = 0; nrow < nsigpairs; nrow++) 
            {
                CompareInfo.sigClusterPair scp = (CompareInfo.sigClusterPair) sigClusterPairs.get(nrow);
                
                if (scp.cluster1.size() >= nmaxcol)
	        {
                    nmaxcol = scp.cluster1.size();
	        }
	        
	        // make rhs profiles, and add each row to rhss
	        ArrayList rhs = new ArrayList(); // rhs of one row, each element is a profile record
	        for (int clusterIndex=0; clusterIndex < scp.rightClusters.size(); clusterIndex++)
	        {
	            ArrayList rightCluster = (ArrayList) (scp.rightClusters.get(clusterIndex));
	            for (int pindex=0; pindex< rightCluster.size();pindex++)
	            {
	                STEM_DataSet.ProfileRec rightProfile = (STEM_DataSet.ProfileRec) ((STEM_DataSet.ProfileRec)rightCluster.get(pindex)).clone();
	                rhs.add(rightProfile);
	            }
	        }
	        rhss.add(rhs); 
	        
	        if (rhs.size() >= nmaxcol)
	        {
                    nmaxcol = rhs.size();
	        }
            }

            CompareInfo.sigClusterPair[] sigrowsArray = new CompareInfo.sigClusterPair[nsigpairs];
            
            for (int nrowindex = 0; nrowindex < sigrowsArray.length; nrowindex++)
            {
                sigrowsArray[nrowindex] = (CompareInfo.sigClusterPair) sigClusterPairs.get(nrowindex);
            }
            
            
            
            //skip sorting
            
            double dheight =20;
            if (nsigpairs > 0)
            {
                dheight = Math.max(Math.min((SCREENHEIGHT-BUFFER)*2/nsigpairs,(SCREENWIDTH*2/3-BUFFER)/(nmaxcol+1)),20);
            }
            double clusterGeneNum = 0;
            double rightClusterGeneNum = 0;
            for (int nrow = 0; nrow < nsigpairs; nrow++) // row (outer) loop
            {
                double dxoffset, dyoffset;

                final CompareInfo.sigClusterPair scp = (CompareInfo.sigClusterPair) sigrowsArray[nrow];
                
                
                
                clusterGeneNum = 0;
                int ncols1 = scp.cluster1.size();
                //resort columns
                STEM_DataSet.ProfileRec[] clusterProfiles1 = new STEM_DataSet.ProfileRec[ncols1];
                for (int ncolindex = 0; ncolindex < clusterProfiles1.length; ncolindex++)
	        {
	            STEM_DataSet.ProfileRec profile = (STEM_DataSet.ProfileRec) (scp.cluster1).get(ncolindex);
	            clusterProfiles1[ncolindex] = (STEM_DataSet.ProfileRec) profile.clone();
	            clusterGeneNum += leftset.countassignments[profile.nprofileindex];
	        }

                // skip sorting
	        double h1=0;
	        double s1=0;
	        
	        for (int ncolindex = 0; ncolindex < ncols1; ncolindex++) // col loop1
	        {
	           final STEM_DataSet.ProfileRec profile = (STEM_DataSet.ProfileRec) clusterProfiles1[ncolindex];

                   PNode colnode = (PNode) leftnodes[profile.nprofileindex].clone();
                   
	           colnode.scale((dheight)/(colnode.getScale()*colnode.getHeight()));
	           
	           // offsets
	           dxoffset = CAPOFFSET/colnode.getScale();
	           

                   dyoffset = colnode.getHeight()*nrow+CAPOFFSET/colnode.getScale();
	           h1 = colnode.getHeight();
                   s1 = colnode.getScale();
	           double h2 = colnode.getHeight();
                   colnode.translate(dxoffset+colnode.getHeight()*(ncolindex+1), dyoffset);
                   
                   if (ncolindex == 0 && nrow==0)
                   {
                       textXOffset1 = dxoffset + colnode.getHeight()*(ncolindex+1);
                       scale1 = colnode.getScale();
                   }
                   

                   NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
                   nf.setMinimumFractionDigits(2);
                   nf.setMaximumFractionDigits(2);
                   String szcorrLabel = nf.format(scp.dcorrval);
                   //if (!theCompareInfo.ColFlag)
                   //{                                   temporarily remove correlation text
                       szcorrLabel = "";
                   //}
                   PText corrtext = new PText(szcorrLabel);
                   corrtext.setFont(new Font("times",Font.PLAIN,(int) (Math.ceil(colnode.getHeight()/6))));
                   corrtext.translate(4*colnode.getWidth()/7,-1);
                   colnode.addChild(corrtext);

	           final String szIntersect = MAINGUI2.doubleToSz(scp.dmatch[ncolindex])+" of the "+
                                       MAINGUI2.doubleToSz(clusterGeneNum)  
                                        +" genes assigned to Profile "
                                      +profile.nprofileindex+" in the " + szother +
                         " experiment were also assigned to this profile (p-value ="
                         +MAINGUI2.doubleToSz(scp.dpvals[ncolindex])+")";

                   
      
                   ArrayList namelist = (ArrayList) ((ArrayList)scp.inames).clone();
                   final HashSet names = (HashSet) namelist.get(ncolindex);
                   colnode.addInputEventListener(new PBasicInputEventHandler() 
                   {
	              public void mousePressed(PInputEvent event) 
                      {
                         if (event.getButton() == MouseEvent.BUTTON1)
                         {          
                            javax.swing.SwingUtilities.invokeLater(new Runnable() 
                            {
                               public void run() 
                               {                
                                  final ProfileGui pg;
                                  
		     	           pg = new ProfileGui(leftset ,profile.nprofileindex,null,null,-1,null,null,null,leftplotpanel,cf); 

			           pg.addWindowListener(new WindowAdapter()
			           {
			             public void windowClosing(WindowEvent we)
			             {
			                pg.dispose();
			             }

			     	     public void windowOpened(WindowEvent we)
				     {
			                pg.repaint();
			             }
			          });

                                  pg.setLocation(20,50);        
                                  //pg.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);  //removed from 1.3.11
                                  pg.setSize(new Dimension(SCREENWIDTH,SCREENHEIGHT));
			          pg.setVisible(true);   
			       }
		           });  
		         }
	              }
	           });    
                   canvas.getLayer().addChild(colnode);       
	        } // end of col loop1                
                
                
                ArrayList rhs = (ArrayList) (rhss.get(nrow));
                
                int ncols2 = rhs.size();
                //resort columns 2
                STEM_DataSet.ProfileRec[] clusterProfiles2 = new STEM_DataSet.ProfileRec[ncols2];
                for (int ncolindex = 0; ncolindex < clusterProfiles2.length; ncolindex++)
	        {
	            STEM_DataSet.ProfileRec profile = (STEM_DataSet.ProfileRec)rhs.get(ncolindex);
	            clusterProfiles2[ncolindex] = (STEM_DataSet.ProfileRec) profile.clone();
	            
	        }

                // skip sorting
	        
	        // RHS
	        // for comparison clusters and profiles (as clusters of 1 profile)
	        for (int ncolindex = 0; ncolindex < ncols2; ncolindex++)
	        {   
	        
	            // count number of genes in this cluster
	            rightClusterGeneNum = 0;
	            CompareInfo.sigClusterPair scpp = (CompareInfo.sigClusterPair) sigClusterPairs.get(nrow);
	            
	            int profileFlag = 0; //0:cluster, 1:profile
	            
	            // find the correct cluster (rightCluster) of the profile being plotted
	            // so that genes in the cluster can be counted later
	            int cnt = ncolindex+1;
	            ArrayList rightCluster = (ArrayList) (scpp.rightClusters.get(0));
	            for (int clusterindex=0; clusterindex < scpp.rightClusters.size();clusterindex++)
	            {   
	                ArrayList rcluster = (ArrayList) (scpp.rightClusters.get(clusterindex));
	                if (rcluster.size() < cnt)
	                {
	                    cnt -= rcluster.size();
	                }
	                else
	                {
	                    rightCluster = (ArrayList) (scpp.rightClusters.get(clusterindex));
	                    break;
	                }
	            }
	            
	            if (rightCluster.size()==1)
	            {
	                profileFlag = 1;
	            }
	            else
	            {
	                profileFlag = 0;
	            }
	            
	            // count the number of genes in the cluster
	            for (int pindex=0; pindex< rightCluster.size();pindex++)
	            {   
	                STEM_DataSet.ProfileRec rightProfile = (STEM_DataSet.ProfileRec)((STEM_DataSet.ProfileRec)rightCluster.get(pindex)).clone();
	                rightClusterGeneNum += rightset.countassignments[rightProfile.nprofileindex];
	            }
	        
	           
	           
	           final STEM_DataSet.ProfileRec profile = (STEM_DataSet.ProfileRec) clusterProfiles2[ncolindex];
                   PNode colnode = (PNode) rightnodes[profile.nprofileindex].clone();
                   

                   dxoffset =1/colnode.getScale()*SCREENWIDTH/5*2+2*CAPOFFSET/colnode.getScale();
                   
                   
	           //dyoffset = h1*nrow+CAPOFFSET/s1;
                   double h2 = colnode.getHeight();
	           colnode.scale((dheight)/(colnode.getScale()*colnode.getHeight()));
	           
	           dyoffset = colnode.getHeight()*nrow+CAPOFFSET/colnode.getScale();
                   colnode.translate(dxoffset+colnode.getHeight()*(ncolindex+1), dyoffset);
                   
                   if (ncolindex == 0 && nrow==0)
                   {
                       textXOffset2 =  dxoffset + colnode.getHeight()*(ncolindex+1);
                       scale2 = colnode.getScale();
                   }                   
                   
                   
                  
                   
                   String szLabel = ((int) scp.dmatch[ncolindex])+";"+MAINGUI2.doubleToSz(scp.dpvals[ncolindex]);
                   PText text = new PText(szLabel);
                   text.setFont(new Font("times",Font.PLAIN,(int) (Math.ceil(colnode.getHeight()/6))));
                   text.translate(1,3*colnode.getHeight()/4);
                   colnode.addChild(text);

                   NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
                   nf.setMinimumFractionDigits(2);
                   nf.setMaximumFractionDigits(2);
                   String szcorrLabel = nf.format(scp.dcorrval);
                   //if (!theCompareInfo.ColFlag)     temporarily remove correlation text
                   //{
                       szcorrLabel = "";
                   //}
                   PText corrtext = new PText(szcorrLabel);
                   corrtext.setFont(new Font("times",Font.PLAIN,(int) (Math.ceil(colnode.getHeight()/6))));
                   corrtext.translate(4*colnode.getWidth()/7,-1);
                   colnode.addChild(corrtext);

                   String clusterOrProfile = "";
                   if (profileFlag == 0) // cluster
                   {
                       clusterOrProfile = "This cluster has "+rightClusterGeneNum+"genes, "+MAINGUI2.doubleToSz(scp.dmatch[ncolindex])+" of them are also assigned to the compared cluster in the "+ szother + " experiment (p-value ="
                         +MAINGUI2.doubleToSz(scp.dpvals[ncolindex])+")";
                   }
                   else
                   {
                       clusterOrProfile = MAINGUI2.doubleToSz(scp.dmatch[ncolindex]) + " genes in this profile are also assigned to the compared cluster in the" + szother + " experiment (p-value ="
                         +MAINGUI2.doubleToSz(scp.dpvals[ncolindex])+")";
                   }

	           final String szIntersect = clusterOrProfile;

                   
                   ArrayList namelist = (ArrayList) ((ArrayList)scp.inames).clone();
                   final HashSet names = (HashSet) namelist.get(ncolindex);
                   colnode.addInputEventListener(new PBasicInputEventHandler() 
                   {
	              public void mousePressed(PInputEvent event) 
                      {
                         if (event.getButton() == MouseEvent.BUTTON1)
                         {          
                            javax.swing.SwingUtilities.invokeLater(new Runnable() 
                            {
                               public void run() 
                               {                
                                  final ProfileGui pg;
                                  
		     	           pg = new ProfileGui(rightset ,profile.nprofileindex,null,
						 names,profile.nprofileindex,szIntersect,null,null,rightplotpanel,cf); 
                                  
			           pg.addWindowListener(new WindowAdapter()
			           {
			             public void windowClosing(WindowEvent we)
			             {
			                pg.dispose();
			             }

			     	     public void windowOpened(WindowEvent we)
				     {
			                pg.repaint();
			             }
			          });

                                  pg.setLocation(20,50);        
                                  //pg.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);  //removed from 1.3.11
                                  pg.setSize(new Dimension(SCREENWIDTH,SCREENHEIGHT));
			          pg.setVisible(true);   
			       }
		           });  
		         }
	              }
	           });    
                   canvas.getLayer().addChild(colnode);       
	        } // end of col loop 2
	        
            } // end of row loop
           
           originaltext.scale(scale1);
           comparetext.scale(scale2);
           originaltext.translate(textXOffset1,0);
           comparetext.translate(textXOffset2,0);
           canvas.getLayer().addChild(comparetext);
           canvas.getLayer().addChild(originaltext);
           

           final CompareGui thisFrame = this;
           thisFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);//exit the app when close this window
            
        } //end of compare clusters condition
        
        
        else // compare profiles condition
        {
            canvas = getCanvas();
            canvas.getLayer().removeAllChildren();
       
            int nmaxcol = 0;
            float frecwidth = (float) 150.0; 
            int OFFSET = 0;
            ArrayList sigrows;
            final STEM_DataSet rowset, colset;
            PNode[] rownodes, colnodes;

            String sztexttop, sztextleft; 
            synchronized (swaplock)
            {
                if (bswap)
                {
                    sigrows = theCompareInfo.sigrowsswap;
                    colset = theCompareInfo.origset;
                    rowset = theCompareInfo.comparesetfmnel;
	            rowplotpanel = theCompareInfo.compareframe.thegeneplotpanel;
	            colplotpanel = thegeneplotpanel;
                    colnodes = profilenodes;
	            rownodes = theCompareInfo.compareframe.profilenodes;
                    sztexttop = ST.szDataFileDEF.split("_")[0];
            	    sztextleft = ST.szCompareDEF.split("_")[0];
                    szother = ST.szCompareDEF.split("_")[0];
                }
                else
                {
                    sigrows = theCompareInfo.sigrows;
                    rowset = theCompareInfo.origset;
                    colset = theCompareInfo.comparesetfmnel;
                    rownodes = profilenodes;
	            colnodes = theCompareInfo.compareframe.profilenodes;
                    colplotpanel = theCompareInfo.compareframe.thegeneplotpanel;
	            rowplotpanel = thegeneplotpanel;
                    sztextleft = ST.szDataFileDEF.split("_")[0];
                    sztexttop = ST.szCompareDEF.split("_")[0];
                    szother = ST.szDataFileDEF.split("_")[0];
                }
            }

            PText comparetext = new PText(sztexttop);
            comparetext.setFont(new Font("times",Font.PLAIN,14));
            comparetext.translate(3.0*SCREENWIDTH/18.0,0);
            canvas.getLayer().addChild(comparetext);

            PText originaltext = new PText(sztextleft);
            originaltext.setFont(new Font("times",Font.PLAIN,14));
            originaltext.translate(0,SCREENHEIGHT/2.0);
            originaltext.rotate(-Math.PI/2);
            canvas.getLayer().addChild(originaltext);

            int nsigrows = sigrows.size();
            if (nsigrows >= 2)
            {
                PText comparetext2 = new PText(sztexttop);
                comparetext2.setFont(new Font("times",Font.PLAIN,14));
                comparetext2.translate(12.0*SCREENWIDTH/18.0,0);
                canvas.getLayer().addChild(comparetext2);

                PText originaltext2 = new PText(sztextleft);
                originaltext2.setFont(new Font("times",Font.PLAIN,14));
                originaltext2.translate(SCREENWIDTH/2.0+CAPOFFSET,SCREENHEIGHT/2.0);
                originaltext2.rotate(-Math.PI/2);
                canvas.getLayer().addChild(originaltext2);
            }

            for (int nrow = 0; nrow < nsigrows; nrow++)
            {
                CompareInfo.CompareInfoRow cir = (CompareInfo.CompareInfoRow) sigrows.get(nrow);         
                if (cir.sigprofiles.size() >= nmaxcol)
	        {
                    nmaxcol = cir.sigprofiles.size();
	        }
            }

            CompareInfo.CompareInfoRow[] sigrowsArray = new CompareInfo.CompareInfoRow[nsigrows];
            for (int nrowindex = 0; nrowindex < sigrowsArray.length; nrowindex++)
            {
                sigrowsArray[nrowindex] = (CompareInfo.CompareInfoRow) sigrows.get(nrowindex);
            }
           
            if (nsort == SORTSIG)
            {
                Arrays.sort(sigrowsArray, new SigRowComparator()); 
            }
            else if (nsort == SORTSUPRISE)
            {
                Arrays.sort(sigrowsArray, new SupriseRowComparator()); 
            }

            double dheight =20;
            if (nsigrows > 0)
            {
                dheight = Math.max(Math.min((SCREENHEIGHT-BUFFER)*2/nsigrows,(SCREENWIDTH/2-BUFFER)/(nmaxcol+1)),20);
            }

            for (int nrow = 0; nrow < nsigrows; nrow++)
            {
                double dxoffset, dyoffset;

                final CompareInfo.CompareInfoRow cir = (CompareInfo.CompareInfoRow) sigrowsArray[nrow];
                PNode node = (PNode) rownodes[cir.nprofile].clone();




	        node.scale((dheight)/(node.getScale()*node.getHeight()));
                if (nrow < Math.ceil(nsigrows/2.0))
	        {
                    dxoffset = CAPOFFSET/node.getScale();
                    dyoffset = node.getHeight()*nrow+CAPOFFSET/node.getScale();
                }
                else
	        {
                    dxoffset =1/node.getScale()*SCREENWIDTH/2+2*CAPOFFSET/node.getScale();
                    dyoffset = node.getHeight()*(nrow-Math.ceil(nsigrows/2.0))+CAPOFFSET/node.getScale();
	        }
            
                String szcountLabel = "" +(int) rowset.countassignments[cir.nprofile];
                PText counttext = new PText(szcountLabel);
                counttext.setFont(new Font("times",Font.PLAIN,(int) (Math.ceil(node.getHeight()/6))));
                counttext.translate(1,node.getHeight()-node.getHeight()/4);
                node.addChild(counttext);
                
                double h1 = node.getHeight();
                double w1 = node.getWidth();
                node.translate(dxoffset, dyoffset);
                if (nrow == 0)
	        {
	            double dwidth = (float)(node.getWidth()*node.getScale());
	            double dtheight = Math.ceil(nsigrows/2.0)*(float)(node.getHeight()*node.getScale());
                    PNode line = PPath.createRectangle((float) (CAPOFFSET+dwidth+1.5*node.getScale()),
                                                 (float) CAPOFFSET,
                                                 (float) (6*node.getScale()), 
						   (float) (dtheight));                   
                    line.setPaint(ST.buttonColor);
                    canvas.getLayer().addChild(line);
	        }
                else if (nrow == Math.ceil(nsigrows/2.0))
	        {
	            double dscale = node.getScale();
	            double dwidth = (float)(node.getWidth()*dscale);

	            double dtheight = (nsigrows-Math.ceil(nsigrows/2.0))*(float)(node.getHeight()*dscale);
                    PNode linemid = PPath.createRectangle((float) (2*CAPOFFSET+dwidth+SCREENWIDTH/2.0+1.5*dscale),
						    (float) CAPOFFSET,(float) (6*dscale),
                                                    (float) (dtheight));                   
                    linemid.setPaint(ST.buttonColor);
                    canvas.getLayer().addChild(linemid);

                    float flast  =(float)((CAPOFFSET/dscale+10+node.getHeight()*(nmaxcol+1))*dscale);
	            float fmid = (float)(flast + (SCREENWIDTH/2.0+CAPOFFSET/dscale-flast)/2);

                    PNode line = PPath.createRectangle(fmid, 
				 (float) 5,(float) (6*dscale), (float) (SCREENHEIGHT-BUFFER+15));                   
                    line.setPaint(ST.lightBlue);
                    canvas.getLayer().addChild(line);
	        }

                node.addInputEventListener(new PBasicInputEventHandler() 
                {
	             public void mousePressed(PInputEvent event) 
                    {
                        if (event.getButton() == MouseEvent.BUTTON1)
                        {                         
                            final ProfileGui pg;
		             pg = new ProfileGui(rowset,cir.nprofile,null,null,-1,null,null,null,rowplotpanel,cf);
   		             pg.addWindowListener(new WindowAdapter()
		             {
		                 public void windowClosing(WindowEvent we)
		                 {
			              pg.dispose();
		                 }

		                 public void windowOpened(WindowEvent we)
		                 {
			              pg.repaint();
		                 }
		             });

		             pg.setLocation(20,50);        
		             //pg.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);  //removed for 1.3.11
		             pg.setSize(new Dimension(SCREENWIDTH,SCREENHEIGHT));
		             pg.setVisible(true);           
		        }    
	            }
	        });      

                canvas.getLayer().addChild(node);

                ArrayList compareprofiles = cir.sigprofiles;

                int ncols = compareprofiles.size();
                //resort columns
                CompareInfo.CompareInfoRec[] cirecArray = new CompareInfo.CompareInfoRec[ncols];
                for (int ncolindex = 0; ncolindex < cirecArray.length; ncolindex++)
	        {
	            cirecArray[ncolindex] = (CompareInfo.CompareInfoRec) compareprofiles.get(ncolindex);
	        }

                if (nsort == SORTSIG)
	        {
                    Arrays.sort(cirecArray, new SigComparator()); 
	        }  
                else if (nsort == SORTSUPRISE)
	        {
                    Arrays.sort(cirecArray, new SupriseComparator(cir.nprofile)); 
	        }

                for (int ncolindex = 0; ncolindex < ncols; ncolindex++)
	        {
	           final CompareInfo.CompareInfoRec cirec = (CompareInfo.CompareInfoRec) cirecArray[ncolindex];
                   PNode colnode = (PNode) colnodes[cirec.nprofile].clone();
	           double h2 = colnode.getHeight();

                   
                   //find the position for this profile node
                   //colnode.setGlobalScale(node.getGlobalScale());
                   colnode.scale((dheight)/(colnode.getScale()*colnode.getHeight()));
                   
                   if (nrow < Math.ceil(nsigrows/2.0))
	           {
                       dxoffset = CAPOFFSET/colnode.getScale();
                       dyoffset = colnode.getHeight()*nrow+CAPOFFSET/colnode.getScale();
                   }
                   else
	           {
                       dxoffset =1/colnode.getScale()*SCREENWIDTH/2+2*CAPOFFSET/colnode.getScale();
                       dyoffset = colnode.getHeight()*(nrow-Math.ceil(nsigrows/2.0))+CAPOFFSET/colnode.getScale();
	           }
                   colnode.translate(10+dxoffset+colnode.getHeight()*(ncolindex+1), dyoffset);
                   
                   
                   
                   
                   String szLabel = ((int) cirec.dmatch)+";"+MAINGUI2.doubleToSz(cirec.dpval);
                   PText text = new PText(szLabel);
                   text.setFont(new Font("times",Font.PLAIN,(int) (Math.ceil(colnode.getHeight()/6))));
                   text.translate(1,3*colnode.getHeight()/4);
                   colnode.addChild(text);

                   NumberFormat nf = NumberFormat.getInstance(Locale.ENGLISH);
                   nf.setMinimumFractionDigits(2);
                   nf.setMaximumFractionDigits(2);
                   String szcorrLabel = nf.format(cirec.dcorrval);
                   if (!theCompareInfo.ColFlag)
                   {
                       szcorrLabel = "";
                   }
                   PText corrtext = new PText(szcorrLabel);
                   corrtext.setFont(new Font("times",Font.PLAIN,(int) (Math.ceil(colnode.getHeight()/6))));
                   corrtext.translate(4*colnode.getWidth()/7,-1);
                   colnode.addChild(corrtext);

	           final String szIntersect = MAINGUI2.doubleToSz(cirec.dmatch)+" of the "+
                                       MAINGUI2.doubleToSz(rowset.countassignments[cir.nprofile])  
                                        +" genes assigned to Profile "
                                      +cir.nprofile+" in the " + szother +
                         " experiment were also assigned to this profile (p-value ="
                         +MAINGUI2.doubleToSz(cirec.dpval)+")";

                   colnode.addInputEventListener(new PBasicInputEventHandler() 
                   {
	              public void mousePressed(PInputEvent event) 
                      {
                         if (event.getButton() == MouseEvent.BUTTON1)
                         {          
                            javax.swing.SwingUtilities.invokeLater(new Runnable() 
                            {
                               public void run() 
                               {                
                                  final ProfileGui pg;
		     	           pg = new ProfileGui(colset ,cirec.nprofile,null,
						 cirec.inames,cir.nprofile,szIntersect,null,null,colplotpanel,cf); 

			           pg.addWindowListener(new WindowAdapter()
			           {
			             public void windowClosing(WindowEvent we)
			             {
			                pg.dispose();
			             }

			     	     public void windowOpened(WindowEvent we)
				     {
			                pg.repaint();
			             }
			          });

                                  pg.setLocation(20,50);        
                                  //pg.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);  //removed from 1.3.11
                                  pg.setSize(new Dimension(SCREENWIDTH,SCREENHEIGHT));
			          pg.setVisible(true);   
			       }
		           });  
		         }
	              }
	           });    
                   canvas.getLayer().addChild(colnode);       
	        }
            }

           PNode supriseButton = PPath.createRectangle((float) 0.0,(float) 0.0,frecwidth,(float) 18.0);
           supriseButton.translate(4*(SCREENWIDTH+OFFSET)/5-100,SCREENHEIGHT-65);
           PText thesupriseText = new PText("Order By Correlation");
           thesupriseText.setFont(new Font("times",Font.PLAIN,12));
           thesupriseText.translate(23,2);
           supriseButton.setPaint(ST.buttonColor);
           supriseButton.addChild(thesupriseText);   
           canvas.getLayer().addChild(supriseButton);
           supriseButton.addInputEventListener(new PBasicInputEventHandler() 
           {
               public void mousePressed(PInputEvent event) 
               {
	           nsort = SORTSUPRISE;
	           drawcomparemain();
                   //repaint();
	       }
           });


           PImage helpButton = new PImage(Util.getImageURL("Help24.gif"));
           canvas.getLayer().addChild(helpButton);
           helpButton.translate(SCREENWIDTH-70,SCREENHEIGHT-68);
           final CompareGui thisFrame = this;
           thisFrame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);//exit the app when close this window
           helpButton.addInputEventListener(new PBasicInputEventHandler() 
           {
               public void mousePressed(PInputEvent event) 
               {
                   javax.swing.SwingUtilities.invokeLater(new Runnable() 
                   {
                       public void run() 
                       {
                           JDialog helpDialog = new JDialog(thisFrame, "Help", false);
                           Container theHelpDialogPane = helpDialog.getContentPane();
        
                           helpDialog.setBackground(Color.white);
                           theHelpDialogPane.setBackground(Color.white);
		            String szMessage;

                           if (rowset.bkmeans)
		            {
		                szMessage = 
                                "The display shows to the left of the yellow line a K-means cluster from one of the data "+
                                "sets.  To the right of the yellow line are K-means clusters from the other data set for which a significant "+
                                "number of genes assigned to the cluster on left were also assigned to it.  "+
                                "The p-value of the number of genes in the intersection is computed using the hypergeometric "+
                                "distribution based on the number of total genes assigned to each of the two clusters, and the "+
	                        "total number of genes on the array.\n\n"+
                                "The display is split in half by a blue line, there is no difference between clusters to the left and right "+
                                "of the blue line.  The clusters in the display can be reordered based on the significance through the "+
                                "'Order by Significance' option.  Within "+
                                "each row the clusters are ordered in decreasing order of significance.  The rows are reordered based on "+
                                "decreasing significance of the most significant intersection for the row.  'Order by Correlation' is "+
                                "similar but instead of re-ordering based on significance, the re-ordering is done based on correlation, "+
                                "so that one can quickly identify dissimilar pairs of clusters.  'Order by ID' is the default method to "+
                                "order clusters, based on increasing cluster ID.  "+
                                "'Swap Rows and Columns' swaps which data set has clusters to the left of the yellow line and "+
                                "which data set has clusters organized in columns to the right of the yellow line.\n\n"+
                                "Note also the main interface is zoomable and pannable, "+
                                "hold the right button down to zoom or the left to pan while moving the mouse.";
		          }
		          else
		          {
                                szMessage = "The display shows to the left of the yellow line a profile from one of the data "+
                                "sets.  To the right of the yellow line are profiles from the other data set for which a significant "+
                                "number of genes assigned to the profile on left were also assigned to it.  "+
                                "The p-value of the number of genes in the intersection is computed using the hypergeometric "+
                                "distribution based on the number of total genes assigned to each of the two profiles, and the "+
	                        "total number of genes on the array.\n\n"+
                                "The display is split in half by a blue line, there is no difference between profiles to the left and right "+
                                "of the blue line.  The profiles in the display can be reordered based on the significance through the "+
                                "'Order by Significance' option.  Within "+
                                "each row the profiles are ordered in decreasing order of significance.  The rows are reordered based on "+
                                "decreasing significance of the most significant intersection for the row.  'Order by Correlation' is "+
                                "similar but instead of re-ordering based on significance, the re-ordering is done based on correlation, "+
                                "so that one can quickly identify dissimilar pairs of profiles.  'Order by ID' is the default method to "+
                                "order profiles, based on increasing profile ID.  "+
                                "'Swap Rows and Columns' swaps which data set has profiles to the left of the yellow line and "+
                                "which data set has profiles organized in columns to the right of the yellow line.\n\n"+
                                "Note also the main interface is zoomable and pannable, "+
                                "hold the right button down to zoom or the left to pan while moving the mouse.";
		          }
	
                                    
                          JTextArea textArea = new JTextArea(szMessage,9,60);
                          textArea.setLineWrap(true);
                          textArea.setWrapStyleWord(true);
                          textArea.setBackground(Color.white);
                          textArea.setEditable(false);

		          ImageIcon ii;

		          if (rowset.bkmeans)
		          {
		            ii = Util.createImageIcon("p6_2.png");
		          }
		          else
		          {
                            ii = Util.createImageIcon("p38_13.png");
		          }

                          JLabel jl = new JLabel(ii);
                          theHelpDialogPane.add(jl);

		          JPanel psl = new JPanel();
		          psl.setLayout(new SpringLayout());
		          psl.setBackground(Color.white);
		          psl.add(jl);
		          JScrollPane jsp2 = new JScrollPane(textArea);
		          psl.add(jsp2);
                          SpringUtilities.makeCompactGrid(psl,2,1,2,2,2,2);
		          JScrollPane jsp =new JScrollPane(psl);

		          theHelpDialogPane.add(jsp);
                          theHelpDialogPane.setSize(800,600);
                          theHelpDialogPane.validate();
                          helpDialog.setLocation(thisFrame.getX()+25,thisFrame.getY()+25);
       
                          helpDialog.setSize(800,600);                  
                          helpDialog.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
                          
                          helpDialog.setVisible(true);
		        }
	              });
	         }
              });          

              PNode idButton = PPath.createRectangle((float) 0.0,(float) 0.0,frecwidth,(float) 18.0);
              idButton.translate(2*(SCREENWIDTH+OFFSET)/5-100,SCREENHEIGHT-65);
              PText theidText = new PText("Order By Profile ID");
              theidText.setFont(new Font("times",Font.PLAIN,12));
              theidText.translate(25,2);
              idButton.setPaint(ST.buttonColor);
              idButton.addChild(theidText);   
              idButton.addInputEventListener(new PBasicInputEventHandler() 
              {
                public void mousePressed(PInputEvent event) 
                {
	            nsort = SORTID;
	            drawcomparemain();
       	   }
              });
              canvas.getLayer().addChild(idButton);       
              canvas.getLayer().addChild(idButton);

              PNode sigButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) frecwidth,(float) 18.0);
              sigButton.translate(3*(SCREENWIDTH+OFFSET)/5-100,SCREENHEIGHT-65);
              PText thesigText = new PText("Order By Significance");
              thesigText.setFont(new Font("times",Font.PLAIN,12));
              thesigText.translate(20,2);
              sigButton.setPaint(ST.buttonColor);
              sigButton.addChild(thesigText);   
              canvas.getLayer().addChild(sigButton);
              sigButton.addInputEventListener(new PBasicInputEventHandler() 
              {
                public void mousePressed(PInputEvent event) 
                {
	            nsort = SORTSIG;
	            drawcomparemain();

	          }
              });
    
             PNode swapButton = PPath.createRectangle((float) 0.0,(float) 0.0,frecwidth,(float) 18.0);
             swapButton.translate(1*(SCREENWIDTH+OFFSET)/5-100,SCREENHEIGHT-65);
             PText theswapText = new PText("Swap Rows and Columns");
             theswapText.setFont(new Font("times",Font.PLAIN,12));
             theswapText.translate(6,2);
             swapButton.setPaint(ST.buttonColor);
             swapButton.addChild(theswapText);   
             canvas.getLayer().addChild(swapButton);

             swapButton.addInputEventListener(new PBasicInputEventHandler() 
             {
                public void mousePressed(PInputEvent event) 
                {
	            synchronized(swaplock)
	            {
	                bswap = !bswap;
	            }

	            drawcomparemain();
	        }
             });
        } // end of compare profiles condition

    }
}
