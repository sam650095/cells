let editdata = [];
document.addEventListener("DOMContentLoaded", async function () {
  const csrftoken = getCookie("csrftoken");
  const preload_results = await fetchAPI(
    "/api/preloadidentifythegates",
    0,
    csrftoken
  );
  console.log(preload_results);
  showdropdown(preload_results.adata_list);
});
function showdropdown(adata_list) {
  const preloadul = document.getElementById("preloadul");

  preloadul.innerHTML = "";
  for (let i = 0; i < adata_list.length; i++) {
    const li = document.createElement("li");
    li.innerHTML = `
                  <a onclick="changeselect('${adata_list[i]}')"
                  class="block px-4 py-2 hover:bg-gray-100 dark:hover:bg-gray-600 dark:hover:text-white m-2"
              >${adata_list[i]}</a>
          `;
    preloadul.appendChild(li);
  }
}
async function changeselect(newText) {
  document.getElementById("selected").textContent = newText;
  const csrftoken = getCookie("csrftoken");
  const form = document.getElementById("identifythegatesform");
  const formData = new FormData(form);
  formData.append(
    "chosen_adata",
    document.getElementById("selected").textContent
  );
  const choseadata = await fetchAPI("/api/choseadata", formData, csrftoken);
  document.getElementById("chosen_results").textContent = choseadata.result;
  document.getElementById("processbutton").classList.remove("hidden");
}
async function processbtn(event) {
  event.preventDefault();
  const csrftoken = getCookie("csrftoken");
  toggleLoading(true, "processbutton");
  const identifythegates_result = await fetchAPI(
    "/api/identifythegates",
    0,
    csrftoken
  );
  console.log(identifythegates_result);
  createTable(identifythegates_result.gate_df);
  toggleLoading(false, "processbutton");

  document.getElementById("btnbox").classList.remove("hidden");
  document.getElementById("clearbtn").classList.remove("hidden");
  document.getElementById("confirmbtn").classList.remove("hidden");
  document.getElementById("nextbtn").classList.remove("hidden");
}
function createTable(data) {
  const table = document.getElementById("nametable");
  const thead = table.querySelector("thead tr");
  const tbody = document.getElementById("nametbody");

  thead.innerHTML = '<th scope="col" class="px-6 py-3">#</th>';
  tbody.innerHTML = "";

  Object.keys(data[0]).forEach((key) => {
    const th = document.createElement("th");
    th.textContent = key;
    th.scope = "col";
    th.className = "px-6 py-3";
    thead.appendChild(th);
  });
  console.log(data);
  editedData = JSON.parse(JSON.stringify(data));

  data.forEach((row, index) => {
    const tr = document.createElement("tr");
    tr.className = "bg-white border-b hover:bg-gray-50 cursor-pointer";

    const tdIndex = document.createElement("td");
    tdIndex.textContent = index + 1;
    tdIndex.className = "px-6 py-4";
    tr.appendChild(tdIndex);

    Object.values(row).forEach((value) => {
      const td = document.createElement("td");
      td.textContent = value !== null ? value : "N/A";
      td.className = "w-1/3 px-6 py-4";
      td.onclick = () => handleClick(td, index);
      tr.appendChild(td);
    });

    tbody.appendChild(tr);
  });
}
function handleClick(td, index) {
  if (td.querySelector("input")) {
    return;
  }
  const div = document.createElement("div");
  div.className = "flex justify-between";
  const currentValue = td.textContent;
  td.textContent = "";
  const input = document.createElement("input");
  input.type = "text";
  input.value = currentValue;
  input.className =
    "my-2 w-36 text-sm text-gray-900 border border-gray-300 rounded-lg bg-gray-50 focus:ring-sky-900 focus:border-sky-900";

  div.appendChild(input);
  const btn = document.createElement("button");
  btn.textContent = "confirm";
  btn.className =
    "rounded-md bg-sky-700 py-2 px-4 text-white text-sm font-semibold my-3";

  div.appendChild(btn);
  td.appendChild(div);
  input.focus();
}
