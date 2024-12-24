#include "ems.hpp"
#include "ems1.hpp"
#include "ems2.hpp"
#include "recursEMS.hpp"
#include "parEMS.hpp"
#include <getopt.h>
#include <cstdlib>

using namespace std;

void usage(const char* argv)
{
    std::cout << "Usage: " << argv << " [OPTIONS] <input-sequence-file>" << std::endl;
    std::cout << "\t-s <version>  Possible versions: 1, 2, 4, 5" << std::endl;
    std::cout << "\t-l <l>        Length (l) of (l,d) motif" << std::endl;
    std::cout << "\t-d <d>        Maximum edit distance (d) of (l,d) motif" << std::endl;
    std::cout << "\t-t <int>      Number of threads" << std::endl;
    exit(-1);
}

int stoi_s(const std::string& s, const int base = 10)
{
    std::size_t pos{};
    try
    {
        const int i = std::stoi(s, &pos, base);
        // std::cout << "std::stoi(" << std::quoted(s) << "): "
        //     << i << "; pos: " << pos << '\n';
        return i;
    }
    catch (std::invalid_argument const& ex)
    {
        std::cerr << "std::invalid_argument::what(): " << ex.what() << '\n';
    }
    catch (std::out_of_range const& ex)
    {
        std::cout << "std::out_of_range::what(): " << ex.what() << '\n';
        const long long ll{std::stoll(s, &pos)};
        std::cout << "std::stoll(" << std::quoted(s) << "): " << ll
            << "; pos: " << pos << '\n';
    }
    catch (...)
    {
        std::cerr << "stoi_s: unknown exception\n";
    }
    throw std::runtime_error{"stoi_s: should not reach here"};
}

int main(int argc, char** argv)
{
    std::cout << "Start the Execution" << std::endl;

    std::string version = "5";
    std::string input;
    Params params{};
    params.num_threads = omp_get_max_threads();
    params.l = 1;
    params.d = 1;

    int option_char;
    while ((option_char = getopt(argc, argv, "l:d:s:t:")) != -1)
    {
        switch (option_char)
        {
        case 'l': params.l = stoi_s(optarg);
            break;
        case 'd': params.d = stoi_s(optarg);
            break;
        case 't': params.num_threads = stoi_s(optarg);
            break;
        case 's': version = string(optarg);
            break;
        case '?': usage(argv[0]);
            break;
        default: ;
        }
    }

    if (params.l <= params.d)
    {
        std::cerr << "l must be greater than d" << std::endl;
        exit(-1);
    }

    if (params.num_threads < 1)
    {
        std::cerr << "Number of threads must be greater than 0. Setting to MAX" << std::endl;
        params.num_threads = omp_get_max_threads();
    }


    if (argc - optind < 1)
    {
        usage(argv[0]);
    }
    input = string(argv[optind]);
    // version = "5"; // we force this as parMES

    if (version == "1")
    {
        cout << "Executing EMS1 " << endl;
        Ems1 ems(input, params.l, params.d, params);
        ems.searchWriteMotifs(params);
    }
    else if (version == "2")
    {
        cout << "Executing EMS2 " << endl;
        Ems2 ems2(input, params.l, params.d, params);
        ems2.searchWriteMotifs(params);
    }
    else if (version == "4")
    {
        cout << "Executing recursEMS " << endl;
        recursEms recems(input, params.l, params.d, params);
        recems.searchWriteMotifs(params);
    }
    else if (version == "5")
    {
        cout << "Executing ParEMS " << endl;
        ParEMS parems(input, params.l, params.d, params);
        parems.searchWriteMotifs(params);
    }
    else
    {
        std::cout << "Unknown version: " << version << std::endl;
        exit(-1);
    }
}
