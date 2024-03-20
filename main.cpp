#include "main.h"

using namespace std;

int read_dataset(vector<simple_sparse_vector> *data_vector, vector<int> *label_vector)
{
    uint dimension = 0;
    string filename = "mnist1.scale";
    // string filename = "data.dat";

    // Start a timer
    // long startTime = get_runtime();

    // OPEN DATA FILE
    // =========================
    std::ifstream data_file(filename.c_str());

    if (!data_file.good())
    {
        std::cerr << "error w/ " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    // Read SVM-Light data file
    // ========================
    int num_examples = 0;
    std::string buf;
    while (getline(data_file, buf))
    {
        // ignore lines which begin with #
        if (buf[0] == '#')
            continue;
        // Erase what comes after #
        size_t pos = buf.find('#');
        if (pos < buf.size())
        {
            buf.erase(pos);
        }
        // replace ':' with white space
        int n = 0;
        for (size_t pos = 0; pos < buf.size(); ++pos)
            if (buf[pos] == ':')
            {
                n++;
                buf[pos] = ' ';
            }
        // read from the string
        std::istringstream is(buf);
        int label = 0;
        is >> label;
        if (label != 1)
        {
            label = -1;
            // std::cerr << "Error reading SVM-light format. Abort." << std::endl;
            // exit(EXIT_FAILURE);
        }
        // printf("number of some %d \n", num_examples);
        label_vector->push_back(label);
        simple_sparse_vector instance(is, n);
        data_vector->push_back(instance);
        num_examples++;
        uint cur_max_ind = instance.max_index() + 1;
        if (cur_max_ind > dimension)
            dimension = cur_max_ind;
    }

    data_file.close();
    printf("read data done!........\n");
    printf("number of sample: %d \n", num_examples);
    return num_examples;
}

int main()
{
    vector<simple_sparse_vector> data_vector = vector<simple_sparse_vector>();
    vector<int> label_vector = vector<int>();
    int dim = read_dataset(&data_vector, &label_vector);
    vector<vector<double> > K(dim, vector<double>(dim, 0.0));
    int f_feat = 780;

    // for (int i =0;i<dim;i++)
    // {
    //     for (int j =0;j<dim;j++)
    //     {
    //         vector<double> r(f_feat, 0.0);
    //         data_vector.at(i).put(&r);
    //         vector<double> q(f_feat, 0.0);
    //         data_vector.at(j).put(&q);
    //         double f = _MEB_kernel(r, label_vector.at(i), q, label_vector.at(j), i, j, 1, 1);
    //         K.at(i)[j] = f;
    //     }
    // }
    // _save_K_file("Map.txt", &K);

    _load_K_file("Map.txt", &K);

    vector<double> z;

    _cacl_z_vector(&K, &z);

    bpcg *algo = new bpcg(dim, K, z);
    bpcg_result res = algo->optimize(20);
    vector<double> alpha = res.alpha;
    vector<int> I_active = res.I_active;

    printf("-------- SVM ---------\n");
    svm *svc = new svm(&alpha, &I_active);
    svc->eval(&data_vector, &label_vector);
    // delete svc;

    delete algo;
    delete &res;

    // delete &K;
    // delete &z;

    return (EXIT_SUCCESS);
}
