import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.*;
import java.util.stream.IntStream;

public class AntColony {

    static ThreadMXBean bean = ManagementFactory.getThreadMXBean();

    /* define constants */
    private double c = 1.0;
    private double alpha = 1;
    private double beta = 5;
    private double evaporation = 0.5;
    private double Q = 500;
    private double antFactor = 0.8;
    private double randomFactor = 0.01;

    private int maxIterations = 1000;

    private int numberOfCities;
    private int numberOfAnts;
    private int[][] graph;
    private double trails[][];
    private List<Ant> ants = new ArrayList<>();
    private static Random random = new Random();
    private double probabilities[];

    private int currentIndex;

    private int[] bestTourOrder;
    private double bestTourLength;
    static double probs[] = null;

    static int numberOfTrials = 15;
    static int MAXINPUTSIZE = 10;
    static String ResultsFolderPath = "/home/nicolocker/Results/"; // pathname to results folder
    static FileWriter resultsFile;
    static PrintWriter resultsWriter;

    public static void main(String[] args) {

        verifyWorks();
        System.out.println("\n");

        System.out.println("Running first full experiment...");
        runFullExperiment("AntColony1-CircularCost.txt");     //change all 3 to RandomCost, EuclideanCost or CircularCost depending on which is being used
        System.out.println("Running second full experiment...");
        runFullExperiment("AntColony2-CircularCost.txt");
        System.out.println("Running third full experiment...");
        runFullExperiment("AntColony3-CircularCost.txt");
    }

    public static void verifyWorks() {
        AntColony.GenerateRandomCostMatrix(10);
        AntColony.GenerateRandomEuclideanCostMatrix(10, 20);
        AntColony.GenerateRandomCircularGraphCostMatrix(10, 40);
    }


    public static void runFullExperiment(String resultsFileName) {

        try {
            resultsFile = new FileWriter(ResultsFolderPath + resultsFileName);
            resultsWriter = new PrintWriter(resultsFile);
        } catch (Exception e) {
            System.out.println("*****!!!!!  Had a problem opening the results file " + ResultsFolderPath + resultsFileName);
            return; // not very foolproof... but we do expect to be able to create/open the file...
        }

        ThreadCpuStopWatch TrialStopwatch = new ThreadCpuStopWatch(); // for timing an individual trial
        double prevTime = 0;

        resultsWriter.println("#NumOfVertices        AvgTime         Best Tour Cost"); // # marks a comment in gnuplot data
        resultsWriter.flush();

        for (int inputSize = 2; inputSize <= MAXINPUTSIZE; inputSize++) {
            System.out.println("Running test for input size " + inputSize + " ... ");
            System.gc();

            long batchElapsedTime = 0;

            int[][] matrix;
            for (long trial = 0; trial < numberOfTrials; trial++) {
                //matrix = AntColony.GenerateRandomCostMatrix(inputSize);      //** uncomment to test this one, comment the below
                //matrix = AntColony.GenerateRandomEuclideanCostMatrix(inputSize, 50);   //** uncomment to test this one, comment the above
                matrix = AntColony.GenerateRandomCircularGraphCostMatrix(inputSize, 50); //** uncomment to test this one, comment the above
                TrialStopwatch.start(); // *** uncomment this line if timing trials individually
                AntColony antColony = new AntColony(matrix);
                antColony.startAntOptimization();
                batchElapsedTime += TrialStopwatch.elapsedTime(); // *** uncomment this line if timing trials individually
            }

            double averageTimePerTrialInBatch = (double) batchElapsedTime / (double) numberOfTrials;

            double ratio = 0;
            if (prevTime > 0) {
                ratio = averageTimePerTrialInBatch / prevTime;
            }

            prevTime = averageTimePerTrialInBatch;

            /* print data for this size of input */
            resultsWriter.printf("%6d  %20.2f %13.2f\n", inputSize, averageTimePerTrialInBatch, ratio); // might as well make the columns look nice
            resultsWriter.flush();
            System.out.println(" ....done.");

        }
    }

    public static int[][] GenerateRandomCostMatrix ( int edgeCost){
        int[][] matrix = new int[edgeCost][edgeCost];

        for (int i = 0; i < edgeCost; i++) {
            //  System.out.println("++++");

            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                    //   System.out.println(matrix[i][j]);

                } else {
                    int temp = random.nextInt(100);
                    matrix[i][j] = temp;
                    //    System.out.println(matrix[i][j]);
                }
            }
        }
        System.out.println("Generate Random Cost Matrix:");
        printMatrix(matrix);

        return matrix;
    }

    public static int[][] GenerateRandomEuclideanCostMatrix ( int numOfVerticies, int coordinate){
        Vertex[] vertice = new Vertex[numOfVerticies];

        for (int i = 0; i < numOfVerticies; i++) {
            vertice[i] = new Vertex(random.nextInt(coordinate - 1), random.nextInt(coordinate - 1), i);
        }

        int[][] matrix = new int[numOfVerticies][numOfVerticies];

        for (int i = 0; i < numOfVerticies; i++) {
            //System.out.println("++++");
            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                    //System.out.println(matrix[i][j]);
                } else {
                    matrix[i][j] = vertice[i].distance(vertice[j]);
                    matrix[j][i] = vertice[j].distance(vertice[i]);
                    //System.out.println(matrix[i][j]);
                }
            }
        }

        System.out.println("Generate Random Euclidean Cost Matrix:");
        printMatrix(matrix);

        return matrix;
    }

    public static int[][] GenerateRandomCircularGraphCostMatrix ( int numOfVerticies, int radius){
        Vertex[] vertice = new Vertex[numOfVerticies];

        double angle = 360 / numOfVerticies;
        double curAngle = 0;
        ArrayList<Vertex> sortVertices = new ArrayList<Vertex>();

        for (int i = 0; i < numOfVerticies; i++) {
            double radian = curAngle * Math.PI / 180;
            double x = Math.cos(radian) * radius;
            double y = Math.sin(radian) * radius;
            sortVertices.add(new Vertex(x, y, i));
            curAngle = curAngle + angle;
        }

        Collections.shuffle(sortVertices);
        vertice = sortVertices.toArray(vertice);

        int[][] matrix = new int[numOfVerticies][numOfVerticies];
        for (int i = 0; i < numOfVerticies; i++) {
            //System.out.println("++++");
            for (int j = 0; j <= i; j++) {
                if (i == j) {
                    matrix[i][j] = 0;
                    //System.out.println(matrix[i][j]);
                } else {
                    matrix[i][j] = vertice[i].distance(vertice[j]);
                    matrix[j][i] = vertice[j].distance(vertice[i]);
                    //System.out.println(matrix[i][j]);
                }
            }
        }

        System.out.println("Generate Circular Graph Cost Matrix:");
        printMatrix(matrix);

        return matrix;
    }

    public static void printMatrix ( int[][] matrix){
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[i].length; j++) {
                System.out.printf("%5d", matrix[i][j]);
            }
            System.out.println();
        }
        System.out.println("\n");
    }

    /*https://github.com/eugenp/tutorials/tree/master/algorithms-genetic/src/main/java/com/baeldung/algorithms/ga/ant_colony*/
    public AntColony(int[][] matrix) {
        graph = matrix;
        //graph = generateRandomMatrix(noOfCities);
        numberOfCities = graph.length;
        numberOfAnts = (int) (numberOfCities * antFactor);

        trails = new double[numberOfCities][numberOfCities];
        probabilities = new double[numberOfCities];
        IntStream.range(0, numberOfAnts)
                .forEach(i -> ants.add(new Ant(numberOfCities)));
    }

    public void startAntOptimization() {
        IntStream.rangeClosed(1, 3)
                .forEach(i -> {
                    System.out.println("Attempt #" + i);
                    solve();
                });
    }

    public int[] solve() {
        setupAnts();
        clearTrails();
        IntStream.range(0, maxIterations)
                .forEach(i -> {
                    moveAnts();
                    updateTrails();
                    updateBest();
                });
        System.out.println("Best tour length: " + (bestTourLength - numberOfCities));
        System.out.println("Best tour order: " + Arrays.toString(bestTourOrder));
        return bestTourOrder.clone();
    }

    private void setupAnts() {
        IntStream.range(0, numberOfAnts)
                .forEach(i -> {
                    ants.forEach(ant -> {
                        ant.clear();
                        ant.visitCity(-1, random.nextInt(numberOfCities));
                    });
                });
        currentIndex = 0;
    }

    private void moveAnts() {
        IntStream.range(currentIndex, numberOfCities - 1)
                .forEach(i -> {
                    ants.forEach(ant -> ant.visitCity(currentIndex, selectNextCity(ant)));
                    currentIndex++;
                });
    }

    private int selectNextCity(Ant ant) {
        int t = random.nextInt(numberOfCities - currentIndex);
        if (random.nextDouble() < randomFactor) {
            OptionalInt cityIndex = IntStream.range(0, numberOfCities)
                    .filter(i -> i == t && !ant.visited(i))
                    .findFirst();
            if (cityIndex.isPresent()) {
                return cityIndex.getAsInt();
            }
        }
        calculateProbabilities(ant);
        double r = random.nextDouble();
        double total = 0;
        for (int i = 0; i < numberOfCities; i++) {
            total += probabilities[i];
            if (total >= r) {
                return i;
            }
        }

        throw new RuntimeException("There are no other cities");
    }

    public void calculateProbabilities(Ant ant) {
        int i = ant.trail[currentIndex];
        double pheromone = 0.0;
        for (int l = 0; l < numberOfCities; l++) {
            if (!ant.visited(l)) {
                pheromone += Math.pow(trails[i][l], alpha) * Math.pow(1.0 / graph[i][l], beta);
            }
        }
        for (int j = 0; j < numberOfCities; j++) {
            if (ant.visited(j)) {
                probabilities[j] = 0.0;
            } else {
                double numerator = Math.pow(trails[i][j], alpha) * Math.pow(1.0 / graph[i][j], beta);
                probabilities[j] = numerator / pheromone;
            }
        }
    }

    private void updateTrails() {
        for (int i = 0; i < numberOfCities; i++) {
            for (int j = 0; j < numberOfCities; j++) {
                trails[i][j] *= evaporation;
            }
        }
        for (Ant a : ants) {
            double contribution = Q / a.trailLength(graph);
            for (int i = 0; i < numberOfCities - 1; i++) {
                trails[a.trail[i]][a.trail[i + 1]] += contribution;
            }
            trails[a.trail[numberOfCities - 1]][a.trail[0]] += contribution;
        }
    }

    private void updateBest() {
        if (bestTourOrder == null) {
            bestTourOrder = ants.get(0).trail;
            bestTourLength = ants.get(0)
                    .trailLength(graph);
        }
        for (Ant a : ants) {
            if (a.trailLength(graph) < bestTourLength) {
                bestTourLength = a.trailLength(graph);
                bestTourOrder = a.trail.clone();
            }
        }
    }

    private void clearTrails() {
        IntStream.range(0, numberOfCities)
                .forEach(i -> {
                    IntStream.range(0, numberOfCities)
                            .forEach(j -> trails[i][j] = c);
                });
    }
}

