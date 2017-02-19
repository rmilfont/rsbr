import java.io.*;
import java.util.*;

import com.sun.corba.se.impl.orbutil.graph.*;
import util.*;
import rbr.*;
import paths_gen.*;
import sbr.*;
import util.Graph;

public class rsbr {

    public static void main(String[] args) {

        String topologyFile = args[0];
        String volumePath = null;
        boolean shouldMerge = true;
        int dim = 5;
        int rows = dim;
        int columns = dim;
        double[] faultyPercentages = {0.0};

        if (args.length == 4) {
            rows = Integer.parseInt(args[0]);
            columns = Integer.parseInt(args[1]);
        }

        for (double faultyPercentage : faultyPercentages) {

            if (!hasEnoughEdges(rows, columns, faultyPercentage))
                break;

            System.out.println("Generating graph");
            Graph graph = (topologyFile != null) ?
                    FromFileGraphBuilder.generateGraph(topologyFile) :
                    RandomFaultyGraphBuilder.generateGraph(rows, columns, (int)Math.ceil(faultyPercentage*numberOfEdges(rows, columns)));

            System.out.println("Isolated?: " + graph.hasIsolatedCores());

            System.out.println(" - SR Section");
            GraphRestrictions restrictions = SBRSection(graph);

            System.out.println("Paths Computation");
            ArrayList<ArrayList<Path>> allMinimalPaths = new PathFinder(graph, restrictions).pathsComputation();
            System.out.println(allMinimalPaths.size());

            System.out.println(" - Paths Selection Section");
            RBR rbr = new RBR(graph);

            Map<Path, Double> volumes = null;
            if (volumePath != null) {
                File commvol = new File(volumePath);
                if (commvol.exists()) {
                    System.out
                            .println("Getting volumes from " + volumePath);
                    volumes = communicationVolume(allMinimalPaths, commvol, graph);
                }
            }

            ArrayList<ArrayList<Path>> chosenPaths = selectPaths(allMinimalPaths, graph, volumes);

            Vertex src = new Vertex("1.0");
            Vertex dst = new Vertex("3.3");
            ArrayList<Path> wantedBag = wantedBag(src,dst, chosenPaths);
            wantedBag.clear();

            List<Path> allCriticalPaths = new PathFinder(graph, restrictions).allPathsBetweenVertices(src, dst);
            Collections.sort(allCriticalPaths, new ByLatency(graph));

            wantedBag.add(allCriticalPaths.get(0));

            System.out.println("Vertex src: " + src.name() + " Vertex dst: " + dst.name());
            System.out.println(" - Quantity of paths: " + allCriticalPaths.size());
            System.out.println("1st path: " + allCriticalPaths.get(0).toString());

            System.out.println(" - RBR Section");
            RBRSection(shouldMerge, graph, chosenPaths, rbr, "custom" + "_" +  "01");
            System.out.println(" -----");


            for(int i = 1; i < allCriticalPaths.size(); i++){
                wantedBag.add(allCriticalPaths.get(i));
                RBRSection(shouldMerge, graph, chosenPaths, rbr, "custom"+ "_" + String.format("%02d",wantedBag.size()));
            }

            for(Path path : allCriticalPaths){
                    System.out.println("Size: " + path.size() + " Latency: " + latencyOf(graph,path) + " Path: " + path.toString());
            }

            int size6 =0, size8 =0, size10 = 0, size12 = 0, size14 = 0, size16 = 0, size18 = 0, size20 = 0;
            for(Path path : allCriticalPaths){
                if(path.size() == 6)
                    size6++;
                if(path.size() == 8)
                    size8++;
                if(path.size() == 10)
                    size10++;
                if(path.size() == 12)
                    size12++;
                if(path.size() == 14)
                    size14++;
                if(path.size() == 16)
                    size16++;
                if(path.size() == 18)
                    size18++;
                if(path.size() == 20)
                    size20++;
            }
            if((size6 + size8 + size10 + size12 + size14 + size16 + size18 + size20) == allCriticalPaths.size())
                System.out.println("Sizes: " + size6 + "(6)  " + size8 + "(8)  " + size10 + "(10)  " + size12 + "(12)  " + size14 + "(14)  " + size16 + "(16)  " + size18 + "(18)" + size20 + "(20)");
            else
                System.out.println("Error!!!");

            Set<Edge> edges = new HashSet<>();

            for(int i = 0; i < allCriticalPaths.size(); i++){
                Set<Edge> currentEdges = NewEdges(graph, allCriticalPaths.get(i), edges);
                //System.out.println(currentEdges.size());
                edges.addAll(currentEdges);
                //System.out.println(i+1 +" path " + pathSize + " tamanho set: " + edges.size() + " peso: " + t);
                System.out.println( getPathsWeight(allCriticalPaths, i)/edges.size());
            }


        }
    }

    private static Set<Edge> NewEdges(Graph graph, Path currentPath, Set<Edge> edges) {
        Set<Edge> currentEdges = edgesFromPath(graph, currentPath);
        Set<Edge> intersection = new HashSet<>(edges);
        intersection.retainAll(currentEdges);
        currentEdges.removeAll(intersection);
        return currentEdges;
    }

    private static double getPathsWeight(List<Path> allCriticalPaths, int i) {
        double pathSize = 0.0;
        for(int k = 0; k < i+1; k++){
            pathSize += allCriticalPaths.get(k).size() - 1;
        }

        pathSize = (pathSize) * (1.0d/(i+1));
        return pathSize;
    }

    private static Set<Edge> edgesFromPath(Graph graph, Path path) {
        Set<Edge> edges = new HashSet<>();
        for(int j = 0; j < path.size()-1; j++){
            Edge edge = graph.adjunct(path.get(j), path.get(j+1));
            edges.add(edge);
        }
        return edges;
    }


    private static List<List<Vertex>> setCriticalPairs(Graph graph, GraphRestrictions restrictions, int qntPairs) {
        List<List<Vertex>> criticalPairs = new ArrayList<>();
        while (qntPairs != 0) {
            boolean flag = false;
            int vertexIndexSrc = (int) (Math.random() * ((double) graph.getVertices().size()));
            int vertexIndexDst = (int) (Math.random() * ((double) graph.getVertices().size()));

            Vertex src = new Vertex(graph.getVertices().get(vertexIndexSrc).name());
            Vertex dst = new Vertex(graph.getVertices().get(vertexIndexDst).name());

            // check if are equal
            if (src.equals(dst))
                continue;

            if(graph.adjunct(src, dst) != null)
                continue;


            if(new PathFinder(graph, restrictions).allPathsBetweenVertices(src, dst).size() < 5)
                continue;

            // put in the list
            for(List<Vertex> pairs : criticalPairs) {
                if(pairs.get(0).equals(src) && pairs.get(1).equals(dst) || pairs.get(0).equals(dst) && pairs.get(1).equals(src)){
                    flag = true;
                    break;
                }
            }

            if(!flag) {
                List<Vertex> currentPair = new ArrayList<>();
                currentPair.add(src);
                currentPair.add(dst);
                criticalPairs.add(currentPair);
                qntPairs--;
            }
        }
        return criticalPairs;
    }

    private static ArrayList<Path> wantedBag(Vertex src, Vertex dst, ArrayList<ArrayList<Path>> paths) {
        for(ArrayList<Path> path : paths) {
            if(path.get(0).src().equals(src) && path.get(0).dst().equals(dst)){
                return path;
            }
        }
        return null;
    }

    private static class ByLatency implements Comparator<Path> {
        Graph graph;

        public ByLatency(Graph graph){
            this.graph = graph;
        }

        @Override
        public int compare(Path p0, Path p1) {
            if(latencyOf(graph, p0) < latencyOf(graph, p1)) return -1;
            if(latencyOf(graph, p0) > latencyOf(graph, p1)) return +1;
            return 0;
        }
    }

    public static Double latencyOf(Graph graph, Path path){
        List<Double> TC_i = new ArrayList<>();
        TC_i.add(1.0);
        for(int i = 0; i < path.size()-1; i++) {
          Edge e = graph.adjunct(path.get(i), path.get(i+1));
          TC_i.add(1/(1-e.weight()));
        }
        TC_i.add(1.0);

        JMMatlab jmMatlab = new JMMatlab();
        return jmMatlab.latencia(32, TC_i);
    }

    private static GraphRestrictions SBRSection(Graph graph) {
        BidimensionalSBRPolicy policy = new BidimensionalSBRPolicy();
        SR sbr = new SR(graph,policy);
        System.out.println("Compute the segments");
        sbr.computeSegments();
        System.out.println("Set the restrictions");
        sbr.setrestrictions();
        return sbr.restrictions();
    }

    private static void RBRSection(boolean shouldMerge, Graph graph, ArrayList<ArrayList<Path>> allMinimalPaths, RBR rbr, String fileSuffix) {
        //System.out.println("Regions Computation");
        rbr.addRoutingOptions(allMinimalPaths);
        rbr.regionsComputation();

        if (shouldMerge) {
            //System.out.println("Doing Merge");
            rbr.merge();
            assert rbr.reachabilityIsOk();
        }

        //System.out.println("Making Tables");
        new RoutingTableGenerator(graph, rbr.regions()).doRoutingTable(fileSuffix);
    }

    private static int numberOfEdges(int rows, int columns) {
        return (columns - 1)*rows + (rows - 1)*columns;
    }

    private static ArrayList<ArrayList<Path>> selectPaths(ArrayList<ArrayList<Path>> paths, Graph g, Map<Path, Double> volumes) {
        ArrayList<ArrayList<Path>> chosenPaths = null;
        LinkWeightTracker lwTracker = new LinkWeightTracker(g, volumes);
        int choice = 3;
        switch (choice) {
            case 0: // No selection (all paths)
                chosenPaths = paths;
                break;
            case 1: // Random selection
                chosenPaths = new RandomPathSelector(paths, lwTracker).selection();
                break;
            case 2: // Minimal weight
                chosenPaths = new ComparativePathSelector(paths, new Path.MinWeightComparator(lwTracker), 10, lwTracker).selection();
                break;
            case 3: // Random selection
                chosenPaths = new MinimalPathSelector(paths, lwTracker, g).selection();
                break;
            case 5: // Maximal weight
                chosenPaths = new ComparativePathSelector(paths, new Path.MaxWeightComparator(lwTracker), 2, lwTracker).selection();
        }
        return chosenPaths;
    }

    private static boolean hasEnoughEdges(int rows, int columns, double percentage) {
        int numberOfFaultyEdges = (int) Math.ceil((double) numberOfEdges(rows, columns) * percentage);
        int numberOfGoodEdges = numberOfEdges(rows, columns) - numberOfFaultyEdges;
        int numberOfVertices = rows * columns;

        return (numberOfGoodEdges >= numberOfVertices - 1);
    }

    private static void printResults(ArrayList<ArrayList<Path>> paths, StatisticalAnalyser statistics) {
        double lwAverage = statistics.averageLinkWeight(paths);
        double lwStdDeviation = statistics.standardDeviationLinkWeight(paths);
        double pwAverage = statistics.averagePathWeight(paths);
        double pwStdDeviation = statistics.standardDeviationPathWeight(paths);
        double pnwAverage = statistics.averagePathNormWeight(paths);
        double pnwStdDeviation = statistics.standardDeviationPathNormWeight(paths);
        double ard = statistics.averageRoutingDistance(paths);
        System.out.println("Peso dos caminhos: "+pwAverage+" ("+pwStdDeviation+")");
        System.out.println("Peso normalizado dos caminhos: "+pnwAverage+" ("+pnwStdDeviation+")");
        System.out.println("Peso dos links: "+lwAverage+" ("+lwStdDeviation+")");
        System.out.println("ARD: " + ard);
    }

    private static Map<Path, Double> communicationVolume(ArrayList<ArrayList<Path>> paths, File commvol, Graph graph) {

        int N = graph.columns()*graph.rows();
        double[][] vol = new double[N][N];
        double maxVol = 0;
        Map<Path, Double> pathsVolume = new HashMap<>();

        try {
            Scanner sc = new Scanner(new FileReader(commvol));

            for(int i = 0; i < N; i++) {
                String[] lines = sc.nextLine().split(" \t");
                for(int j = 0; j < N; j++) {
                    vol[i][j] = Double.valueOf(lines[j]);
                    maxVol = (vol[i][j] > maxVol) ? vol[i][j] : maxVol;
                }
            }
            sc.close();

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        for(ArrayList<Path> alp : paths) {
            int i = graph.indexOf(alp.get(0).src());
            int j = graph.indexOf(alp.get(0).dst());
            double volume = vol[i][j];
            for(Path path : alp) {
                pathsVolume.put(path, volume/maxVol);
            }
        }
        return pathsVolume;
    }
}
