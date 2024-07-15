let new_rename_df;
let samplemethodtextnode;
// renameing
async function grabnames() {
  const csrftoken = getCookie("csrftoken");
  const grabname_result = await fetchAPI("/api/grabnames", 0, csrftoken);
  show_name_table(grabname_result.rename_df);
}
function show_name_table(rename_df) {
  const tableBody = document.getElementById("nametbody");

  // Clear existing rows
  tableBody.innerHTML = "";
  new_rename_df = rename_df;
  console.log(new_rename_df.NewName);
  // Add new rows based on the API result
  for (let i = 0; i < rename_df.CurrentName.length; i++) {
    const row = document.createElement("tr");
    row.className = "bg-white border-b hover:bg-gray-50 ";

    row.innerHTML = `
      <th scope="row" class="px-6 py-4 font-medium text-gray-900 whitespace-nowrap ">
        ${i}
      </th>
      <td class="px-6 py-4 font-medium text-gray-900 whitespace-nowrap w-1/3" id="currentname_${i}">
        ${rename_df.CurrentName[i]}
      </td>
      <td class="px-6 py-4 w-1/3" id="newname_${i}">${rename_df.NewName[i]}</td>
      <td class="px-6 py-4 text-right w-1/6">
        <a onclick="editname(${i})" id="editbtn_${i}" class="font-medium text-sky-700 underline-offset-4 hover:underline hover:cursor-pointer" >Edit</a>
        <a onclick="confirmname(${i})" id="confirmbtn_${i}" class="hidden font-medium text-sky-700 underline-offset-4 hover:underline hover:cursor-pointer" >Confirm</a>
      </td>
    `;

    tableBody.appendChild(row);
  }
}
function editname(index) {
  const newname_td = document.getElementById("newname_" + index);
  const old_name = newname_td.textContent;
  newname_td.innerHTML = `
  <input type="text" id="input_${index}" value="${old_name}" class="p-3 w-36 text-sm text-gray-900 border border-gray-300 rounded-lg bg-gray-50 focus:ring-sky-900 focus:border-sky-900">
  `;
  document.getElementById("editbtn_" + index).classList.add("hidden");
  document.getElementById("confirmbtn_" + index).classList.remove("hidden");
}
function confirmname(index) {
  const input = document.getElementById("input_" + index);
  const newname_td = document.getElementById("newname_" + index);
  const currentname_td = document.getElementById("currentname_" + index);
  newname_td.textContent = input.value;
  currentname_td.textContent = input.value;
  new_rename_df.NewName[index] = input.value;
  document.getElementById("editbtn_" + index).classList.remove("hidden");
  document.getElementById("confirmbtn_" + index).classList.add("hidden");
}
async function reanme_confirmbtn() {
  console.log(new_rename_df);
  const csrftoken = getCookie("csrftoken");
  const jsonData = JSON.stringify(new_rename_df);
  const renameresult = await fetchAPI("/api/rename", jsonData, csrftoken);
  loadImage("cluster_result", "clustering_summary.png", "summary-container");
  loadImage("cluster_result", "clustering_heatmap.png", "heatmap-container");
  loadImage("cluster_result", "clustering_leidens.png", "umap-container");
  loadImage("cluster_result", "clustering_ranking.png", "ranking-container");
}
// subclustering
async function grabclusters() {
  const csrftoken = getCookie("csrftoken");
  const grabclusters_result = await fetchAPI("/api/grabclusters", 0, csrftoken);
  console.log(grabclusters_result);
  show_divideclusters(grabclusters_result.clusters_list);
}
function show_divideclusters(clusters) {
  const dropdownList = document.querySelector("#divide_subcluster_dropdown ul");
  dropdownList.innerHTML = "";
  document.getElementById("_").textContent = clusters[0];
  clusters.forEach((cluster, index) => {
    const li = document.createElement("li");
    li.innerHTML = `
        <div class="flex p-2 rounded hover:bg-gray-100">
          <div class="flex items-center h-5">
            <input
              id="cluster-${index}"
              aria-describedby="cluster-text-${index}"
              type="checkbox"
              value="${cluster}"
              checked
              class="clusters_checkbox w-4 h-4 text-blue-600 bg-gray-100 border-gray-300 rounded focus:ring-blue-500 focus:ring-2"
            />
          </div>
          <div class="ms-2 text-sm w-full">
            <label for="cluster-${index}" class="font-bold text-gray-900">
              <div>${cluster}</div>
            </label>
          </div>
        </div>
      `;
    dropdownList.appendChild(li);
  });
}
async function subcluster_confirmbtn() {
  const csrftoken = getCookie("csrftoken");
  const form = document.getElementById("subclusterform");

  const formData = new FormData(form);

  const selectedCheckboxes = document.querySelectorAll(
    ".clusters_checkbox:checked"
  );
  const selectedClusters = Array.from(selectedCheckboxes).map(
    (checkbox) => checkbox.value
  );

  selectedClusters.forEach((cluster) => {
    formData.append("clusters", cluster);
  });

  const subclusters_result = await fetchAPI(
    "/api/subclusters",
    formData,
    csrftoken
  );
  console.log(subclusters_result);
  loadImage("cluster_result", "clustering_summary.png", "summary-container");
  loadImage("cluster_result", "clustering_heatmap.png", "heatmap-container");
  loadImage("cluster_result", "clustering_leidens.png", "umap-container");
  loadImage("cluster_result", "clustering_ranking.png", "ranking-container");
}
// subset
async function preload_subset() {
  const csrftoken = getCookie("csrftoken");
  const preload_subset_result = await fetchAPI(
    "/api/preloadsubset",
    0,
    csrftoken
  );
  console.log(preload_subset_result);
  // data text
  let available_files_result = document.getElementById(
    "available_files_result"
  );
  available_files_result.innerHTML = "";
  let ul = document.createElement("ul");
  ul.classList.add("space-y-1", "text-gray-500", "list-disc", "list-inside");
  preload_subset_result["available_files_result"].forEach((item) => {
    let li = document.createElement("li");
    li.textContent = item;
    ul.appendChild(li);
  });
  available_files_result.appendChild(ul);
  // dropdown show
  show_dropdown(preload_subset_result["clustering_columns"]);
}
function show_dropdown(clustering_col) {
  console.log(clustering_col);
  const subset_sample = document.getElementById("subset_sample");
  const sampleul = document.getElementById("sampleul");

  sampleul.innerHTML = "";
  subset_sample.removeChild(subset_sample.firstChild);

  samplemethodtextnode = document.createTextNode("None");
  subset_sample.insertBefore(samplemethodtextnode, subset_sample.firstChild);
  for (let i = 0; i < clustering_col.length; i++) {
    const li = document.createElement("li");
    li.innerHTML = `
                <a onclick="samplemethodtextnodefix('${clustering_col[i]}')"
                class="block px-4 py-2 hover:bg-gray-100 m-2"
            >${clustering_col[i]}</a>
        `;
    sampleul.appendChild(li);
  }
}
async function samplemethodtextnodefix(newText) {
  const csrftoken = getCookie("csrftoken");
  samplemethodtextnode.nodeValue = newText;
  const samplemethod = document.getElementById("samplemethod");
  samplemethod.value = newText;
  const formData = new FormData();
  formData.append("sample", newText);
  const grab_cluster_subset_result = await fetchAPI(
    "/api/grabclustersubset",
    formData,
    csrftoken
  );
  console.log(grab_cluster_subset_result);
  show_dropdown_subset(grab_cluster_subset_result.cluster);
}
function show_dropdown_subset(clusters) {
  const dropdownList = document.querySelector("#subset_dropdown ul");
  dropdownList.innerHTML = "";

  clusters.forEach((value, index) => {
    const li = document.createElement("li");
    li.innerHTML = `
        <div class="flex p-2 rounded hover:bg-gray-100">
          <div class="flex items-center h-5">
            <input
              id="subset-${index}"
              aria-describedby="subset-text-${index}"
              type="checkbox"
              value="${value}"
              checked
              class="subset-checkbox w-4 h-4 text-blue-600 bg-gray-100 border-gray-300 rounded focus:ring-blue-500 focus:ring-2"
            />
          </div>
          <div class="ms-2 text-sm w-full">
            <label for="subset-${index}" class="font-bold text-gray-900">
              <div>${value}</div>
            </label>
          </div>
        </div>
      `;
    dropdownList.appendChild(li);
    const checkbox = li.querySelector(`#subset-${index}`);
    checkbox.addEventListener("change", logCheckedValues);
  });
  logCheckedValues();
}
function logCheckedValues() {
  const checkedBoxes = document.querySelectorAll(".subset-checkbox:checked");
  const checkedValues = Array.from(checkedBoxes).map((box) => box.value);
  const selectedBox = document.getElementById("subset_rslt");
  selectedBox.textContent = checkedValues.join(", ");
}
async function subset_confirmbtn() {
  const csrftoken = getCookie("csrftoken");
  const form = document.getElementById("subsetform");

  const formData = new FormData(form);
  formData.append("sample", document.getElementById("samplemethod").value);
  const selectedCheckboxes = document.querySelectorAll(
    ".subset-checkbox:checked"
  );
  const selectedCluster = Array.from(selectedCheckboxes).map(
    (checkbox) => checkbox.value
  );

  selectedCluster.forEach((cluster) => {
    formData.append("clusters", cluster);
  });
  formData.append(
    "naming_cluster",
    document.getElementById("name_cluster").value
  );
  const subset_results = await fetchAPI("/api/subset", formData, csrftoken);
  console.log(subset_results);
}
