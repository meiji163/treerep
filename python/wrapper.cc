#include "pybind_func.h"

PYBIND11_MODULE(trep,m){
	py::class_<Graph>(m, "Graph")
		.def(py::init<>())
		.def("add_edge", &Graph::add_edge)
		.def("adj_list", &Graph::adj_list)
		.def("remove_edge", &Graph::remove_edge)
		.def("retract", &Graph::retract)
		.def("remove_vertex", &Graph::remove_vertex)
		.def("save", &Graph::to_mtx);

	m.def("graph_metric", &py_metric);
	m.def("load_graph", &graph_from_mtx);
	m.def("treerep_graph", &py_graph_treerep, 
			py::arg("G"), py::arg("tol")=0.1);
	m.def("treerep", &py_treerep,
			py::arg("metric"), py::arg("N"), py::arg("tol")=0.1);
}

