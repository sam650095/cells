let new_rename_df;
// renameing
async function grabnames() {
  const csrftoken = getCookie("csrftoken");
  const grabname_result = await fetchAPI(
    "/api/grabphenotypesnames",
    0,
    csrftoken
  );
  console.log(grabname_result);
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
        <td class="px-6 py-4 font-medium text-gray-900 w-1/3" id="currentname_${i}">
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
async function rename_confirmbtn() {
  console.log(new_rename_df);
  const csrftoken = getCookie("csrftoken");
  const jsonData = JSON.stringify(new_rename_df);
  const renameresult = await fetchAPI(
    "/api/phenotyperename",
    jsonData,
    csrftoken
  );
  console.log(renameresult);
  loadImage(
    "phenotype_result",
    "phenotyping_summary.png",
    "p-summary-container"
  );
  loadImage(
    "phenotype_result",
    "phenotyping_heatmap.png",
    "p-heatmap-container"
  );
  loadImage("phenotype_result", "phenotyping_leidens.png", "p-umap-container");
  loadImage(
    "phenotype_result",
    "phenotyping_ranking.png",
    "p-ranking-container"
  );
}
async function grabdropphenotype() {
  const csrftoken = getCookie("csrftoken");
  const grab_drop_phenotype_result = await fetchAPI(
    "/api/grabdropphenotype",
    0,
    csrftoken
  );
  console.log(grab_drop_phenotype_result);
  insert_dreop_list(grab_drop_phenotype_result["drop_list"]);
}
function insert_dreop_list(drop_list) {
  const dropdownList = document.querySelector(
    "#divide_dropphenotype_dropdown ul"
  );
  dropdownList.innerHTML = "";

  drop_list.forEach((marker, index) => {
    const li = document.createElement("li");
    li.innerHTML = `
          <div class="flex p-2 rounded hover:bg-gray-100">
            <div class="flex items-center h-5">
              <input
                id="marker-${index}"
                aria-describedby="marker-text-${index}"
                type="checkbox"
                value="${marker}"
                class="drop-checkbox w-4 h-4 text-blue-600 bg-gray-100 border-gray-300 rounded focus:ring-blue-500 focus:ring-2"
              />
            </div>
            <div class="ms-2 text-sm w-full">
              <label for="marker-${index}" class="font-bold text-gray-900">
                <div>${marker}</div>
              </label>
            </div>
          </div>
        `;
    dropdownList.appendChild(li);
  });
}
async function dropphenotype_confirmbtn() {
  const csrftoken = getCookie("csrftoken");

  const form = document.getElementById("dropphenotypeform");

  const formData = new FormData(form);

  const selectedCheckboxes = document.querySelectorAll(
    ".drop-checkbox:checked"
  );
  console.log(selectedCheckboxes);
  const selecteddrop = Array.from(selectedCheckboxes).map(
    (checkbox) => checkbox.value
  );
  selecteddrop.forEach((drop) => {
    formData.append("drops", drop);
  });
  const dropresult = await fetchAPI("/api/dropphenotype", formData, csrftoken);
  console.log(dropresult);
  loadImage(
    "phenotype_result",
    "phenotyping_summary.png",
    "p-summary-container"
  );
  loadImage(
    "phenotype_result",
    "phenotyping_heatmap.png",
    "p-heatmap-container"
  );
  loadImage("phenotype_result", "phenotyping_leidens.png", "p-umap-container");
  loadImage(
    "phenotype_result",
    "phenotyping_ranking.png",
    "p-ranking-container"
  );
}
