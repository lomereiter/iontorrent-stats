import std.stdio;

public class Logger {
    static void warn(string msg) {
        stderr.writeln("[warning] ", msg);
    }
}
