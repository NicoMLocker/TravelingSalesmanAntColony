public class Vertex {
    public int x;
    public int y;
    public int name;

    public Vertex(double a, double b, int n) {
        x = (int) a;
        y = (int) b;
        name = n;
    }

    public int distance(Vertex point) {
        return (int) Math.sqrt((point.x - this.x) * (point.x - this.x) + (point.y - this.y) * (point.y - this.y));
    }
}