#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>

//TO DO PRELAMANJE U OPSTEM SLUCAJU, BRISANJE, SEARCH, ISPIS
using namespace std;

const int nn = 5; //red stabla br pokazivaca
const int order = 5; //2 * floor((2 * nn - 2) / 3) + 1; OVO JE MAX BR STRINGOVA, SAMO KOREN IMA TOLIKO OSTALI
//SE POPUNJAVAJU DO NN-1

struct Node {

private:

	void init() {
		for (int i = 0; i < order; i++) {
			data[i] = nullptr;
			pointers[i] = nullptr;
		}
		pointers[order] = nullptr;
		currElems = 0;
		father = nullptr;
	}

public:
	string* data[order];
	Node* pointers[order + 1];
	int leaf;
	int root;
	int currElems;
	Node* father;
	Node() {
		init();
	}

	Node(const Node& n);
	Node(Node&& n);
	Node& operator=(const Node& n);
	Node& operator=(Node&& n);
	~Node();
	int maxPointers() {
		if (root)
			return 2 * floor((2 * nn - 2) / 3) + 1;
		else return nn;
	}
	int minPointers() {
		if (root) return 2;
		else return ceil((2 * nn - 1) / 3);
	}

	bool isFull() { return currElems == maxPointers() - 1; }

	friend ostream& operator<<(ostream& os, Node* root) {
		queue<Node*>qu;
		qu.push(root);

		Node* gran = new Node;
		gran->data[0] = new string{ "XXXXX" };

		while (!qu.empty()) {
			Node* tmp = qu.front();
			qu.pop();
			
			if (*tmp->data[0] == "XXXXX")os << endl;
			else {
				os << "[";
				for (int i = 0; i < tmp->currElems; i++) {
					os << *tmp->data[i] << " ";
				}
				os << "]   ";

				for (int i = 0; i < tmp->currElems; i++) {
					if (tmp->pointers[i])qu.push(tmp->pointers[i]);
				}
				qu.push(gran);
			}
		}
		return os;
	}
};

Node* findLeafToInsert(Node* root, string d, int& pos) {
	Node* temp = root, * find = root;
	pos = -1;
	int greater = 0;
	while (temp) {
		find = temp;
		for (int i = 0; i < temp->currElems; i++) {
			if (*temp->data[i] > d) {
				temp = temp->pointers[i], pos = i;
				if (!temp)return find;
			}
			else {
				greater++;
			}
		}
		if (greater == temp->currElems)pos = temp->currElems, temp = temp->pointers[temp->currElems]; //ako je veci od svih
	}
	return find;
}

//a b 

void deleteKeys(Node* nod) {
	for (int i = 0; i < nod->currElems; i++) {
		delete nod->data[i];
		nod->data[i] = nullptr;
	}
}

void nodeSeparating(Node* separate, string d, int pos) {//pos nam pozicija prelomnog cvora
	int n = separate->currElems + (separate->root ? 1 : 2);
	int inserted = 0;

	int firstNum = floor((2 * nn - 2) / 3);
	int secondNum = floor((2 * nn - 1) / 3);
	int thirdNum = floor((2 * nn) / 3);

	vector<string> allData;

	for (int i = 0; i < separate->currElems; i++)allData.push_back(*separate->data[i]);
	allData.push_back(d);

	if (separate->root) {//koren nemamo razdelni oko koga prelamamo
		sort(allData.begin(), allData.end());

		deleteKeys(separate);

		Node* firstNode = new Node;
		for (int i = 0; i < firstNum; i++) {
			firstNode->data[i] = new string{ allData[i] };
			firstNode->currElems++;
		}
		separate->data[0] = new string{ allData[firstNum] };
		separate->pointers[0] = firstNode;
		separate->pointers[0]->leaf = 1;
		separate->pointers[0]->root = 0;
		inserted += firstNum + 1;
		firstNode->father = separate;
		separate->currElems -= firstNode->currElems;

		Node* secondNode = new Node;
		for (int i = 0; i < secondNum; i++) {
			secondNode->data[i] = new string{ allData[i + firstNum + 1] };
			secondNode->currElems++;
		}
		inserted += secondNum;
		secondNode->father = separate;
		separate->pointers[1] = secondNode;
		separate->pointers[1]->leaf = 1;
		separate->pointers[1]->root = 0;
		separate->currElems -= secondNode->currElems;

		if (inserted != n) {// ako se uzme 3 granicnik mora da ima 3 cvora
			separate->data[1] = new string{ allData[inserted] };
			inserted++;

			Node* thirdNode = new Node;
			for (int i = 0; i < thirdNum; i++) {
				thirdNode->data[i] = new string{ allData[i + firstNum + secondNum + 2] };
				thirdNode->currElems++;
			}
			thirdNode->father = separate;
			separate->pointers[2] = thirdNode;
			separate->pointers[2]->leaf = 1;
			separate->pointers[2]->root = 0;
			separate->currElems -= thirdNode->currElems;
			inserted += thirdNum;
		}
	}
	else {//PRELAMANJE CVORA KOJI NIJE KOREN?
		allData.push_back(*separate->data[pos]); // prelomni kljuc iz oca
		//uzima desnog, ako nema desnog onda levog brata koji je takodje popunjen
		Node* Brother = (separate->father->pointers[pos + 1] ? separate->father->pointers[pos + 1] : separate->father->pointers[pos - 1]);
		for (int i = 0; i < Brother->currElems; i++)allData.push_back(*Brother->data[i]);
		Node* fath = separate->father;
		sort(allData.begin(), allData.end());

		deleteKeys(separate);
		deleteKeys(Brother);

		Node* firstNode = separate->father->pointers[pos + 1] ? separate : separate->father->pointers[pos - 1];

		for (int i = 0; i < firstNum; i++) {
			firstNode->data[i] = new string{ allData[i] };
		}
		inserted += firstNum + 1;

		Node* secondNode = separate->father->pointers[pos + 1] ? separate->father->pointers[pos + 1] : separate;
		for (int i = 0; i < secondNum; i++) {
			secondNode->data[i] = new string{ allData[i + firstNum + 1] };
		}
		inserted += secondNum;

		if (inserted != n) {// ako se uzme 3 granicnik mora da ima 3 cvora
			inserted++;

			Node* thirdNode = new Node;
			for (int i = 0; i < thirdNum; i++) {
				thirdNode->data[i] = new string{ allData[i + firstNum + secondNum + 2] };
			}
			inserted += thirdNum;
		}


		//napravimo 3 cvora i saljemo sad u oca
		if (separate->father->currElems + 2 < separate->father->maxPointers() - 1) {
			//ubacujemo ta dva kljuca i prepovezujemo
		}
		else {
			//cuvamo pokazivace iz oca pre pos, oca delimo
		}

	}

}

void nodePouring(Node* pour, Node* brother, int pos) {//pos lokacija prelomnog
	vector<string> allData;
	Node* fath = pour->father;
	for (int i = 0; i < pour->currElems; i++)allData.push_back(*pour->data[i]);
	for (int i = 0; i < brother->currElems; i++)allData.push_back(*brother->data[i]);
	allData.push_back(*pour->father->data[pos]);

	sort(allData.begin(), allData.end());

	int mid = allData.size() / 2;
	delete fath->data[pos];
	fath->data[pos] = new string{ allData[mid] };

	deleteKeys(pour);
	deleteKeys(brother);

	for (int i = 0; i < mid; i++)pour->data[i] = new string{ allData[i] };
	for (int i = mid + 1; i < allData.size(); i++)brother->data[i] = new string{ allData[i] };

}

void insertNode(Node* root, string d) {
	int positionOfSon = 0;

	if (root->data[0] == nullptr) {//prvi string umecemo
		root->data[0] = new string{ d };
		root->currElems += 1;
	}
	else {
		Node* place = findLeafToInsert(root, d, positionOfSon);
		if (place->currElems < place->maxPointers()) {//ima mesta u cvoru
			for (int j = place->currElems; j > positionOfSon; j--)place->data[j] = place->data[j - 1];
			place->data[positionOfSon] = new string{ d };
			place->currElems++;
		}
		else {//prelivanje ako nema bracu prelamanje
			int flag = 0;
			if (place->root) {//kod njega nema prelivanja samo prelamanje
				nodeSeparating(place, d, positionOfSon);

			}
			else {
				Node* help = nullptr;
				//prelivanje ako moze, ako ne prelamamo
				if (!place->father->pointers[positionOfSon + 1]->isFull()) {//da li je u opsegu
					help = place->father->pointers[positionOfSon + 1];
					flag = 1;
				}
				else if (!place->father->pointers[positionOfSon - 1]->isFull()) { //da li je u opsegu
					help = place->father->pointers[positionOfSon - 1];
					flag = 1;
				}
				if (flag) {
					nodePouring(place, help, positionOfSon);
				}
				else {
					//prelamanje
				}
			}

		}
	}
}


int main() {
	Node* root = new Node;
	root->leaf = 0;
	root->root = 1;
	insertNode(root, "a");
	insertNode(root, "c");
	insertNode(root, "b");
	insertNode(root, "d");
	insertNode(root, "e");
	insertNode(root, "jk");
	cout << root;
}