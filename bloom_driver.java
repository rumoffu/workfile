package bf;
import java.util.*;
import java.io.*;
public class bloom_driver{

    public static class read
    {
		 /* read class:
		  * - name read1002
		 * - seq AGCT
		 * - kmer_len length of kmers used
		 * - num_kmers total
		 * - #kmers in +
		 * - #kmers in -
		 * - actual_pos T/F
		 * - called_pos T/F*/
    	
    	//constructor
 	   read( String namey, String seqy, int len )
 	   {
 		   name  	= namey;
 		   seq = seqy;
 		   kmer_len = len;
 	   }
 	   String name;
	   String seq;
	   int kmer_len;
	   int num_kmers = 0;
	   int num_pos_kmers = 0;
	   int num_neg_kmers = 0;
	   boolean actual_pos;
	   boolean called_pos;
	   /*
	   public void set_num_kmers(int num_kmersy) {
	       this.num_kmers = num_kmersy;
	   }
	   public int get_num_kmers() {
	       return this.num_kmers;
	   }
	   public void set_num_pos_kmers(int num_pos_kmersy) {
	       this.num_pos_kmers = num_pos_kmersy;
	   }
	   public int get_num_pos_kmers() {
	       return this.num_pos_kmers;
	   }
	   public void set_num_neg_kmers(int num_neg_kmersy) {
	       this.num_neg_kmers = num_neg_kmersy;
	   }
	   public int get_num_neg_kmers() {
	       return this.num_neg_kmers;
	   }
	   public void set_actual_pos(boolean actual_posy) {
	       this.actual_pos = actual_posy;
	   }
	   public boolean get_actual_pos() {
	       return this.actual_pos;
	   }
	   public void set_called_pos(boolean called_posy) {
	       this.called_pos = called_posy;
	   }
	   public boolean get_called_pos() {
	       return this.called_pos;
	   }
	   */
    }
	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		/* 1) read in + and - enhancers and count 
		 * 2) save first 1/4 for testing,
		 * 3) use 3/4 of enhancers to break into kmers and train + and - bloom filters
		 * 4) run in test data 
		 * 5) report # correct and data
		 * */

        
		boolean consoleprint = false;
		double false_positive_probability = 0.1;
		int kmer_size = 0;
		boolean normalize = false;
		int testcase = 0;
		if(args.length == 3)
		{
			kmer_size = Integer.parseInt(args[0]);
			testcase = Integer.parseInt(args[1]);
			normalize = Boolean.parseBoolean(args[2]);
		}
		else
		{
			System.out.println("usage: kmer_size(int) testcase(1-4) normalize(true/false)");
			System.exit(0);
		}
		PrintWriter toFile = new PrintWriter(new FileWriter("sum.txt"));
		toFile.println("kmer_size\ttestcase\tnormalize\t%+right\t%-right\t%avgright\t%phit\t%nhit\truntime\tpos_diff\tneg_diff");
		//for normalize and non-normal
		for(int norm = 0; norm < 2; norm++)
		{
			
		if(norm == 0) normalize = false;
		else if(norm == 1) normalize = true;
		
		//for testcase
		for(int tc = 1; tc <= 4; tc++)
		{
			testcase = tc;
		//for loop for kmer's
		for(int km = 5; km < 21; km++)
		{
	        // Get and store the current time -- for timing
	        long runstart;
	        runstart = System.currentTimeMillis();
	        
			kmer_size = km;
			//if(km == 21) kmer_size = 50;
		Scanner infile1 = new Scanner( new FileReader( "enh_fb.fa" ) );
		Scanner infile2 = new Scanner( new FileReader( "nullseqsi_200_1.fa" ) );
		String title;
		String seq;
		int in_size = 5000;
		read[] pos_reads = new read[in_size];
		read[] neg_reads = new read[in_size];
		int num_pos_reads = 0;
		int num_neg_reads = 0;
		int total_pos_length = 0;
		int total_neg_length = 0;
		
		
		while(infile1.hasNext())
		{//read in positive naively
			title = infile1.nextLine();
			if(infile1.hasNext())
			{
				seq = infile1.nextLine();
				pos_reads[num_pos_reads++] = new read(title.substring(1), seq, kmer_size);
				total_pos_length += seq.length();
			}

			
		}
		while(infile2.hasNext())
		{//read in negative naively
			title = infile2.nextLine();
			if(infile2.hasNext())
			{
				seq = infile2.nextLine();
				neg_reads[num_neg_reads++] = new read(title.substring(1), seq, kmer_size);
				total_neg_length += seq.length();
			}
			
		}
		
		if(normalize)
		{ //ensure to use lower value
			if( num_neg_reads > num_pos_reads) num_neg_reads = num_pos_reads;
			else num_pos_reads = num_neg_reads; 
		}
		//test on 1/4, train on 3/4
		int[] ptest = new int[num_pos_reads / 4];
		int[] ntest = new int[num_neg_reads / 4];
		int[] ptrain = new int[num_pos_reads - num_pos_reads/4];
		int[] ntrain = new int[num_neg_reads - num_neg_reads/4];
		int pcount = 0; int ncount = 0; int pt = 0; int nt = 0;

		switch(testcase)
		{
			case 1: //test on first 1/4, train on rest
				for(int i = 0; i < num_pos_reads / 4; i++)
					ptest[pt++] = i;
				for(int i = 0; i < num_neg_reads / 4; i++)
					ntest[nt++] = i;
				for(int i = num_pos_reads/4; i < num_pos_reads; i++)
					ptrain[pcount++] = i;
				for(int i = num_neg_reads/4; i < num_neg_reads; i++)
					ntrain[ncount++] = i;
				break;
			case 2://test on second 1/4, train on rest
				for(int i = num_pos_reads / 4; i < (2 * num_pos_reads) / 4; i++)
					ptest[pt++] = i;
				for(int i = num_neg_reads / 4; i < (2 * num_neg_reads) / 4; i++)
					ntest[nt++] = i;
				for(int i = 0; i < num_pos_reads / 4; i++)
					ptrain[pcount++] = i;
				for(int j = (2 * num_pos_reads) / 4; j < num_pos_reads; j++)
					ptrain[pcount++] = j;
				for(int i = 0; i < num_neg_reads / 4; i++)
					ntrain[ncount++] = i;
				for(int j = (2 * num_neg_reads) / 4; j < num_neg_reads; j++)
					ntrain[ncount++] = j;
				break;
			case 3://test on third 1/4, train on rest
				for(int i = (2 * num_pos_reads) / 4; i < (3 * num_pos_reads) / 4; i++)
					ptest[pt++] = i;
				for(int i = (2 * num_neg_reads) / 4; i < (3 * num_neg_reads) / 4; i++)
					ntest[nt++] = i;
				for(int i = 0; i < (2 * num_pos_reads) / 4; i++)
					ptrain[pcount++] = i;
				for(int j = (3 * num_pos_reads) / 4; j < num_pos_reads; j++)
					ptrain[pcount++] = j;
				for(int i = 0; i < (2 * num_neg_reads) / 4; i++)
					ntrain[ncount++] = i;
				for(int j = (3 * num_neg_reads) / 4; j < num_neg_reads; j++)
					ntrain[ncount++] = j;
				break;
			case 4://test on fourth 1/4, train on rest
				for(int i = (3 * num_pos_reads) / 4; i < num_pos_reads -1; i++)
					ptest[pt++] = i;
				for(int i = (3 * num_neg_reads) / 4; i < num_neg_reads -1; i++)
					ntest[nt++] = i;
				for(int i = 0; i < (3 * num_pos_reads) / 4; i++)
					ptrain[pcount++] = i;
				for(int i = 0; i < (3 * num_neg_reads) / 4; i++)
					ntrain[ncount++] = i;
				break;
		}
		//real training length is only 3/4 so calculate
		int use_pos_length = 0; int use_neg_length = 0;
		for(int i = 0; i < ptrain.length; i++) 
		{
			use_pos_length += pos_reads[ptrain[i]].seq.length();
		}

		
		for(int i = 0; i < ntrain.length; i++) 
		{
			use_neg_length += neg_reads[ntrain[i]].seq.length();
		}
		
		BloomFilter<String> bloom_filter_pos = new BloomFilter<String>(false_positive_probability, use_pos_length);
		BloomFilter<String> bloom_filter_neg = new BloomFilter<String>(false_positive_probability, use_neg_length);
		
		String kmer = "";
		//training pos and neg filters
		for(int i = 0; i < ptrain.length; i++) 
		{//for training portion of the reads
			for(int j = 0; j <= pos_reads[ptrain[i]].seq.length() - kmer_size; j++) 
			{//for each character make kmer of kmer_size
				kmer = pos_reads[ptrain[i]].seq.substring(j, j + kmer_size);
				bloom_filter_pos.add(kmer);	
				
			}
			pos_reads[ptrain[i]].num_kmers = pos_reads[ptrain[i]].seq.length() - kmer_size + 1;
			pos_reads[ptrain[i]].actual_pos = true;
		}
		for(int i = 0; i < ntrain.length; i++) 
		{//for training portion of the reads
			for(int j = 0; j <= neg_reads[ntrain[i]].seq.length() - kmer_size; j++) 
			{//for each character make kmer of kmer_size
				kmer = neg_reads[ntrain[i]].seq.substring(j, j + kmer_size);
				bloom_filter_neg.add(kmer);	
			}
			neg_reads[ntrain[i]].num_kmers = neg_reads[ntrain[i]].seq.length() - kmer_size + 1;
			neg_reads[ntrain[i]].actual_pos = false;
		}

		//test on 1/4 of saved data - the test set
		int num_pos_right = 0; //number of positive reads called positive
		int num_neg_right = 0; //number of negative reads called negative
		
		double pphit = 0; double pnhit = 0;
		double nphit = 0; double nnhit = 0;
		int num_pkmers = 0; int num_nkmers = 0;
		
		//positive tests
		for(int i = 0; i < ptest.length; i++) 
		{//for test portion of the reads
			for(int j = 0; j <= pos_reads[ptest[i]].seq.length()-kmer_size; j++) 
			{//for each character make kmer of kmer_size
				kmer = pos_reads[ptest[i]].seq.substring(j, j+kmer_size);
				if(bloom_filter_pos.contains(kmer))
				{
					pos_reads[ptest[i]].num_pos_kmers++;
				}
				if(bloom_filter_neg.contains(kmer))
				{
					pos_reads[ptest[i]].num_neg_kmers++;	
				}
			}
			pos_reads[ptest[i]].num_kmers = pos_reads[ptest[i]].seq.length()-kmer_size + 1;
			pos_reads[ptest[i]].actual_pos = true;
			if( pos_reads[ptest[i]].num_pos_kmers - pos_reads[ptest[i]].num_neg_kmers > 0)
			{
				////ppos_neg_diff += pos_reads[ptest[i]].num_pos_kmers - pos_reads[ptest[i]].num_neg_kmers;
				pos_reads[ptest[i]].called_pos = true;
				num_pos_right += 1;
			}
			else //more negative so call negative
			{
				////ppos_neg_diff += pos_reads[ptest[i]].num_pos_kmers - pos_reads[ptest[i]].num_neg_kmers;
				pos_reads[ptest[i]].called_pos = false;
			}
			pphit += pos_reads[ptest[i]].num_pos_kmers;
			pnhit += pos_reads[ptest[i]].num_neg_kmers;
			num_pkmers += pos_reads[ptest[i]].num_kmers;
		}
		//negative tests
		for(int i = 0 ; i < ntest.length; i++) 
		{//for first 1/4 of the reads
			for(int j = 0; j <= neg_reads[ntest[i]].seq.length()-kmer_size; j++) 
			{//for each character make kmer of kmer_size
				kmer = neg_reads[ntest[i]].seq.substring(j, j+kmer_size);
				if(bloom_filter_pos.contains(kmer))
				{
					neg_reads[ntest[i]].num_pos_kmers++;
				}
				if(bloom_filter_neg.contains(kmer))
				{
					neg_reads[ntest[i]].num_neg_kmers++;
				}
			}
			neg_reads[ntest[i]].num_kmers = neg_reads[ntest[i]].seq.length()-kmer_size + 1;
			neg_reads[ntest[i]].actual_pos = false;
			if( neg_reads[ntest[i]].num_pos_kmers - neg_reads[ntest[i]].num_neg_kmers > 0)
			{
	////			npos_neg_diff += neg_reads[ntest[i]].num_pos_kmers - neg_reads[ntest[i]].num_neg_kmers;
				neg_reads[ntest[i]].called_pos = true;
			}
			else //more negative so call negative
			{
	////			npos_neg_diff += neg_reads[ntest[i]].num_pos_kmers - neg_reads[ntest[i]].num_neg_kmers;
				neg_reads[ntest[i]].called_pos = false;
				num_neg_right += 1;
			}
			nphit += neg_reads[ntest[i]].num_pos_kmers;
			nnhit += neg_reads[ntest[i]].num_neg_kmers;
			num_nkmers += neg_reads[ntest[i]].num_kmers;
		}
		String s = "" + kmer_size + "_" + testcase + "_" + normalize + ".txt";
		double ppos_neg_diff = 0; //absolute difference between # kmers called pos and neg
		double npos_neg_diff = 0; //same but for negative tests
		PrintWriter toFile1 = new PrintWriter(new FileWriter(s));
		for(int i = 0; i < ptest.length; i++)
		{
			ppos_neg_diff += pos_reads[ptest[i]].num_pos_kmers - pos_reads[ptest[i]].num_neg_kmers;
			toFile1.println(i + "\t" + pos_reads[ptest[i]].num_pos_kmers + "\t" + pos_reads[ptest[i]].num_neg_kmers + "\t" + pos_reads[ptest[i]].num_kmers);
		}
		toFile1.println("NEGATIVE");
		for(int i = 0; i < ntest.length; i++)
		{
			npos_neg_diff += neg_reads[ntest[i]].num_pos_kmers - neg_reads[ntest[i]].num_neg_kmers;
			toFile1.println(i + "\t" + neg_reads[ntest[i]].num_pos_kmers + "\t" + neg_reads[ntest[i]].num_neg_kmers + "\t" + neg_reads[ntest[i]].num_kmers);
		}
		toFile1.close();

		double average_ppos_neg_difference = ppos_neg_diff / (ptest.length);
		double average_npos_neg_difference = npos_neg_diff / (ntest.length);
		if(consoleprint)
		{
		System.out.println("kmer_size: " + kmer_size);
		System.out.println("#positive right: " + num_pos_right);
		System.out.println("#positive tests: " + ptest.length);
		System.out.println("% positive right: " + num_pos_right / (double) ptest.length);
		System.out.println("#negative right: " + num_neg_right);
		System.out.println("#negtaive tests: " + ntest.length);
		System.out.println("% negative right: " + num_neg_right / (double) ntest.length);
		////System.out.println("average_ppos_neg_difference: " + average_ppos_neg_difference);
		double average_kmers = (total_pos_length + total_neg_length) / (num_pos_reads + num_neg_reads);
		System.out.println("average_kmers per read: " + average_kmers);
		System.out.println("average % right:"+"\t" + (num_pos_right + num_neg_right ) / (double) (ptest.length + ntest.length));
		}
		String data = "";
		data = data + kmer_size + "\t";
		data = data + testcase + "\t";
		data = data + normalize + "\t";
		data = data + num_pos_right / (double) ptest.length + "\t";
		data = data + num_neg_right / (double) ntest.length + "\t";
		data = data + (num_pos_right + num_neg_right ) / (double) (ptest.length + ntest.length)+ "\t"; //avg right
		data = data + pphit / num_pkmers+ "\t";
		data = data + nphit / num_nkmers+ "\t";
      	// Output search time 
        long elapsed = System.currentTimeMillis() - runstart;
        data = data + elapsed + "\t";
        data = data + average_ppos_neg_difference + "\t";
        data = data + average_npos_neg_difference + "\t";
		toFile.println(data);
 
     	
		}//end kmer loop
		}//end testcase loop
		}//end normalize loop
		toFile.close();
	}

}
